# purpose: create 3-layer network data + analysis
# ------------------------------------------------------------------------
#                           ~ N O T E S ~
# ------------------------------------------------------------------------
#  * Required packages
#    - data.table
#    - openxlsx
#    - igraph
#    - tidygraph
#    - dplyr
#    - tidyr
#    - stringr
#    - ggplot2
#    - GGally
#    - ggpubr
#    - SFdegree (this is a custom package saved locally on HHD's drive)
#
#  * The three layers are:
#    - phylogenetic (phylo)
#    - person-person contact (PP)
#    - intrahousehold (HH)
#
#  * We construct 3 types of network, increasing in complexity
#    1) Phylo-only
#    2) Phylo + PP
#    3) Phylo + PP + HH
#
#  * Create the netwrok using following criteria
#    - new phylogenetic network from Garrick
#    - set patristic distance at 1 month
#      (approx: 160 node, 140 edge)
#    - build social ties using all aggregated
#      (social data is the overall aggregated)
#    - social network should contain nodes that are directly connected
#      to the nodes found in the phylogenetic layer
##########################################################################
rm(list = ls())

path.data1 = ".../Linking MAVEN and Sequenced Data"
path.data2 = ".../Contact Tracing"
path.data3 = ".../Specimen Data - MAVEN"
path.data4 = ".../Patient Address History"
path.data5 = ".../data/RData"
path.data6 = ".../Received From UGA"

path.save = ".../output folder"

today_yyyymmdd = gsub("-", "", as.character(Sys.Date()))

days_in_1month = 30.4167

##########################################################################



# set up -----------
### read data ~ linking files
n8753_linked = openxlsx::read.xlsx(paste(path.data1, "/N8753_linked_to_maven.xlsx", sep = ""))
ryan9800 = data.table::fread(paste(path.data1, "/IDsWLinks.csv", sep = ""))
ryan9800 = ryan9800 |> dplyr::filter(!ID %in% c("HELLO", "TEST-123"))


### read data ~ original contact tracing
pp_data = data.table::fread(paste(path.data2, "/ContactTracing NonCluster Links from202003to202112 - v20220217.csv", sep = ""))


### read data ~ intrahousehold
intrahouse = data.table::fread(paste(path.data4, "/FilterMatches_6_13_wOverlapNLocationID.csv", sep = ""))
intrahouse = intrahouse |>
  dplyr::arrange(EVENT_ID_x, EVENT_ID_y) |>
  dplyr::select(EVENT_ID_x, EVENT_ID_y, OVERLAP_START, OVERLAP_END, LOCATION_ID, POSTAL_CODE)


### testing results
spec_data = data.table::fread(paste(path.data3, "/MAVEN Specimens Sub Data from202003to202202 - v20220406.csv", sep = ""))


### read data ~ phylogenetic network's node set file
phylo_node_set = data.table::fread(paste(path.data6, "/Patristic_Distance_Network_Hou_6moMax_NodeSet.csv", sep = ""))


### read data ~ new phylo network from Garrick ... v2022-12-21
meta.orig = data.table::fread(paste(path.data1, "/1. 2. PamBrown/Texas_cumulative2.csv", sep = ""))
phylo_net_with_ct = data.table::fread(paste(path.data6, "/patristic20221221.csv", sep = "")) |>
  dplyr::distinct() |>
  # get collection date ~ 'source'
  dplyr::left_join(
    meta.orig |> dplyr::select(GISAID_name, collection_date),
    by = c("source" = "GISAID_name")
  ) |>
  # get collection date ~ 'target'
  dplyr::left_join(
    meta.orig |> dplyr::select(GISAID_name, collection_date),
    by = c("target" = "GISAID_name"),
    suffix = c(".src", ".tar")
  )

### remove dup'd edges
phylo_net_with_ct = phylo_net_with_ct |>
  dplyr::rowwise() |>
  dplyr::mutate(temp = paste(min(source, target), "-", max(source, target), sep = "")) |>
  dplyr::ungroup() |>
  dplyr::distinct(temp, .keep_all = TRUE) |>
  dplyr::select(-c(temp))

############


# some functions ------
  ### a function to extract sequence "ID" based on some referencing MAVEN Event ID
    # `link_file` should be Ryan's 9800 file (contain Maven Event ID & sequencing ID)
  get_SeqID_from_EventID = function(link_file, id_refernce, eid_reference) {
    ## some input arg check
    if (!all(c("ID", "Maven_Event_ID") %in% colnames(link_file)))
      stop("`ID` and `Maven_Event_ID` should be column names in `link_file`")


    ## first approach to extract the ID
    id1 = link_file |>
      dplyr::filter(ID %in% id_refernce) |>
      dplyr::group_by(ID) |>
      dplyr::filter(any(Maven_Event_ID %in% eid_reference)) |>
      dplyr::ungroup() |>
      dplyr::pull(ID) |>
      unique()

    ## second approach to extract the ID
    id2 = link_file |>
      dplyr::filter(Maven_Event_ID %in% eid_reference) |>
      dplyr::group_by(ID) |>
      dplyr::filter(any(ID %in% id_refernce)) |>
      dplyr::ungroup() |>
      dplyr::pull(ID) |>
      unique()

    ## comapre if result are the same
    c1 = length(id1) == length(id2)
    c2 = all(id1 %in% id2)

    ## return
    if (!(c1 && c2)) {
      warning("Someting might be not working propoerly...\n")
    }
    out = list(
      check_passed = c1 && c2,
      length = length(id1),
      ids = id1
    )
    return(out)
  }


  ### ............................................................................................................
  ### ....... Given an igraph object, randomly samples `n_clusters` and generate the randomly sampled graph ......
  ### ............................................................................................................
  ### output of this function is NULL unless get_coord is set to `TRUE`
  ###   * if get_coord is set to `TRUE` then the sampled sub-network + x,y coordinates are returned in a list object
  ### input of this function are as the following:
  ###   * `g` is an igraph object for which will be plotted
  ###   * `n_clusters` is a numeric scalar, telling the function how many random clusters to sampled
  ###   * `min_csize` is a numeric scalar, indicating the smallest component size that is to be sampled.
  ###   * `s` is the seed value
  ###   * `get_coords` a logical scalar, if TRUE, the function will return the randomly selected graph + its XY coords on the visualization
  plot_ig.rand_cluster = function(g, n_clusters, min_csize = 2, s = 1, get_coords = FALSE, ...){
    # clustering
    cl = igraph::clusters(g)
    cl_id = which(cl$csize >= min_csize)

    # sub plot
    set.seed(s)
    cl_id.rand = sample(cl_id, size = min(n_clusters, length(cl_id)), replace = FALSE)
    g.sub = igraph::delete.vertices(graph = g, v = !(igraph::get.vertex.attribute(g, name = "name") %in% names(cl$membership[cl$membership %in% cl_id.rand])))

    # plot
    coords = igraph::layout.fruchterman.reingold(g.sub)
    igraph::V(g.sub)$x = coords[,1]
    igraph::V(g.sub)$y = coords[,2]

    plot(x = g.sub, vertex.label = NA)

    if(get_coords){
      return(list(graph = g.sub, coord = coords))
    }else{
      return(NULL)
    }
  }

  
####################



# creating the first layer ~ phylogenetic ------
  ### phylo network subset to 1 month
  phylo_net_with_ct.1month  = phylo_net_with_ct |> dplyr::filter(`patristic distance`*365.25 <= days_in_1month*1)
  cat(
    "## SUMMARY [1] ## ",
    sprintf("\n - `phylo_net_with_ct.1month` has a max patristic distance of: %.0f days", max(phylo_net_with_ct.1month$`patristic distance`*365,25)),
    sprintf("\n -                                     the number of edges is: %.0f edges", nrow(phylo_net_with_ct.1month)),
    sprintf("\n -                              the number of unique nodes is: %.0f nodes", length(unique(c(phylo_net_with_ct.1month$source, phylo_net_with_ct.1month$target)))),
    "\n",
    sep = ""
  )


  ### data frame used to extract edge list and node set
  phylo_net.df = phylo_net_with_ct.1month |>
    dplyr::left_join(
      n8753_linked |> dplyr::select(ID, GISAID_name),
      by = c("source" = "GISAID_name")
    ) |>
    dplyr::left_join(
      n8753_linked |> dplyr::select(ID, GISAID_name),
      by = c("target" = "GISAID_name"),
      suffix = c(".src", ".tar")
    )

  ### edge list + node set
  phylo_net.el = phylo_net.df |> dplyr::select(ID.src, ID.tar) |> dplyr::rename(from = ID.src, to = ID.tar)
  phylo_net.ns = data.frame(name = unique(c(phylo_net.df$ID.src, phylo_net.df$ID.tar)))

  ### igraph
  phylo_net.ig = igraph::graph_from_data_frame(d = phylo_net.el, vertices = phylo_net.ns, directed = FALSE)

  ### add any attributes here
  igraph::V(phylo_net.ig)$set = "phyl"
  igraph::E(phylo_net.ig)$set = "phyl"


#############################################

# ID numbers for the 162 nodes -----
  ### GISAID for the 162 nodes in phylo network who also have contact tracing
  GISAID_for_162.pfx = c(phylo_net_with_ct.1month$source,
                         phylo_net_with_ct.1month$target)
  GISAID_for_162.pfx = unique(GISAID_for_162.pfx)
  GISAID_for_162 = sub("^hCoV-19/", "", GISAID_for_162.pfx)


  ### sequence ID
  ID_162 = n8753_linked |> dplyr::filter(GISAID_name %in% GISAID_for_162.pfx) |> dplyr::pull(ID)


  ### MAVEN Event ID
    # have multiple E_ID for single sample --> not de-duplicated in MAVEN
    # have multiple sample for single E_ID --> same person have multiple samples
  EID_for_162.df = ryan9800 |> dplyr::filter(ID %in% ID_162)
  EID_for_162 = EID_for_162.df |> dplyr::pull(Maven_Event_ID) |> unique()

##################################




# creating the social network layer of network ......... -----
  ### create data frame that we will use to extract edge list and node set
  pp_net.df = pp_data |>
    # subset to those directly connected to phylo nodes
    dplyr::filter(Source_EventID %in% EID_for_162 | Exposed_EventID %in% EID_for_162) |>
    dplyr::select(Source_EventID, Exposed_EventID) |>
    # remove dup'd edges
    dplyr::rowwise() |>
    dplyr::mutate(temp = paste(min(Source_EventID, Exposed_EventID), "-", max(Source_EventID, Exposed_EventID), sep = "")) |>
    dplyr::ungroup() |>
    dplyr::distinct(temp, .keep_all = TRUE) |>
    dplyr::select(-c(temp)) |>
    # change ID to match phylo IDs
    dplyr::left_join(
      EID_for_162.df,
      by = c("Source_EventID" = "Maven_Event_ID")
    ) |>
    dplyr::left_join(
      EID_for_162.df,
      by = c("Exposed_EventID" = "Maven_Event_ID"),
      suffix = c(".source", ".exposed")
    ) |>
    dplyr::mutate(
      Source_EventID = ifelse(!is.na(ID.source), ID.source, Source_EventID),
      Exposed_EventID = ifelse(!is.na(ID.exposed), ID.exposed, Exposed_EventID)
    )

  ### edge list + node set
  pp_net.ns = data.frame(name = unique(c(pp_net.df$Source_EventID, pp_net.df$Exposed_EventID)))
  pp_net.el = pp_net.df |>
    dplyr::select(Source_EventID, Exposed_EventID) |>
    dplyr::rename(from = Source_EventID, to = Exposed_EventID) |>
    dplyr::filter(from != to) |>
    # remove dup'd edges
    dplyr::rowwise() |>
    dplyr::mutate(temp = paste(min(from, to), "-", max(from, to), sep = "")) |>
    dplyr::ungroup() |>
    dplyr::distinct(temp, .keep_all = TRUE) |>
    dplyr::select(-c(temp))

  ### igraph object
  pp_net.ig = igraph::graph_from_data_frame(d = pp_net.el, vertices = pp_net.ns, directed = FALSE)


  ### any network attributes go here
  igraph::V(pp_net.ig)$set = "pp"
  igraph::E(pp_net.ig)$set = "pp"


############################################################

# creating the intra-household layer of the network .... ---------
  ### create data frame that we will use to extract edge list and node set
  hh_net.df = intrahouse |>
    # subset to those directly connected to phylo nodes
    dplyr::select(EVENT_ID_x, EVENT_ID_y) |>
    dplyr::filter(EVENT_ID_x %in% EID_for_162 | EVENT_ID_y %in% EID_for_162) |>
    # remove dup'd edges
    dplyr::rowwise() |>
    dplyr::mutate(temp = paste(min(EVENT_ID_x, EVENT_ID_y), "-", max(EVENT_ID_x, EVENT_ID_y), sep = "")) |>
    dplyr::ungroup() |>
    dplyr::distinct(temp, .keep_all = TRUE) |>
    dplyr::select(-c(temp)) |>
    # change ID to match phylo IDs
    dplyr::left_join(
      EID_for_162.df,
      by = c("EVENT_ID_x" = "Maven_Event_ID")
    ) |>
    dplyr::left_join(
      EID_for_162.df,
      by = c("EVENT_ID_y" = "Maven_Event_ID"),
      suffix = c(".source", ".exposed")
    ) |>
    dplyr::mutate(
      EVENT_ID_x = ifelse(!is.na(ID.source), ID.source, EVENT_ID_x),
      EVENT_ID_y = ifelse(!is.na(ID.exposed), ID.exposed, EVENT_ID_y)
    )

  ### edge list + node set
  hh_net.ns = data.frame(name = unique(c(hh_net.df$EVENT_ID_x, hh_net.df$EVENT_ID_y)))
  hh_net.el = hh_net.df |>
    dplyr::select(EVENT_ID_x, EVENT_ID_y) |>
    dplyr::rename(from = EVENT_ID_x, to = EVENT_ID_y) |>
    dplyr::filter(from != to) |>
    # remove dup'd edges
    dplyr::rowwise() |>
    dplyr::mutate(temp = paste(min(from, to), "-", max(from, to), sep = "")) |>
    dplyr::ungroup() |>
    dplyr::distinct(temp, .keep_all = TRUE) |>
    dplyr::select(-c(temp))

  ### igraph object
  hh_net.ig = igraph::graph_from_data_frame(d = hh_net.el, vertices = hh_net.ns, directed = FALSE)


  ### any network attributes go here
  igraph::V(hh_net.ig)$set = "hh"
  igraph::E(hh_net.ig)$set = "hh"


############################################################



# overlay [1] ~ Phy + PP ... -----
  ### merge node sets
  phylpp_net.ns = rbind(
    phylo_net.ns |> dplyr::mutate(set = "phyl"),
    pp_net.ns |> dplyr::mutate(set = "pp")
  )
  phylpp_net.ns = phylpp_net.ns |>
    dplyr::group_by(name) |>
    dplyr::summarise(set = paste(set, collapse = "")) |>
    dplyr::ungroup() |>
    data.frame()

  ### merge edge list
  phylpp_net.el = rbind(
    phylo_net.el |> dplyr::mutate(set = "phyl"),
    pp_net.el |> dplyr::mutate(set = "pp")
  )
  phylpp_net.el = phylpp_net.el |>
    # create an edge ID
    dplyr::rowwise() |>
    dplyr::mutate(temp = paste(min(from, to), "-", max(from, to), sep = "")) |>
    dplyr::ungroup() |>
    # combine set, if appears in multiple set
    dplyr::group_by(temp) |>
    dplyr::mutate(set = paste(set, collapse = "")) |>
    dplyr::ungroup() |>
    # remove dup'd edges
    dplyr::distinct(temp, set, .keep_all = TRUE) |>
    dplyr::select(-c(temp)) |>
    data.frame()

  ### create igraph
  phylpp_net.ig = igraph::graph_from_data_frame(d = phylpp_net.el, vertices = phylpp_net.ns, directed = FALSE)


################################

# overlay [2] ~ PhyPP + HH . -----
  ### merge node sets
  phylpphh_net.ns = rbind(
    phylpp_net.ns,
    hh_net.ns |> dplyr::mutate(set = "hh")
  )
  phylpphh_net.ns = phylpphh_net.ns |>
    dplyr::group_by(name) |>
    dplyr::summarise(set = paste(set, collapse = "")) |>
    dplyr::ungroup() |>
    data.frame()

  ### merge edge list
  phylpphh_net.el = rbind(
    phylpp_net.el,
    hh_net.el |> dplyr::mutate(set = "hh")
  )
  phylpphh_net.el = phylpphh_net.el |>
    # create an edge ID
    dplyr::rowwise() |>
    dplyr::mutate(temp = paste(min(from, to), "-", max(from, to), sep = "")) |>
    dplyr::ungroup() |>
    # combine set, if appears in multiple set
    dplyr::group_by(temp) |>
    dplyr::mutate(set = paste(set, collapse = "")) |>
    dplyr::ungroup() |>
    # remove dup'd edges
    dplyr::distinct(temp, set, .keep_all = TRUE) |>
    dplyr::select(-c(temp)) |>
    data.frame()

  ### create igraph
  phylpphh_net.ig = igraph::graph_from_data_frame(d = phylpphh_net.el, vertices = phylpphh_net.ns, directed = FALSE)


################################



# (1.0) network summary stat ...................................... -----
  ### generic function to print summary stat
  net_summary = function(g) {

    id_in_net = igraph::V(g)$name
    n_id_in_phylo = sum(id_in_net %in% ID_162)
    n_id_NOT_in_phylo = sum(!id_in_net %in% ID_162)

    cat(
      "\n ", sprintf("               N node: %.0f + %.0f (phylo + non-phylo)", n_id_in_phylo, n_id_NOT_in_phylo),
      "\n ", sprintf("                       = %.0f", length(id_in_net)),
      "\n ", sprintf("N edge (density in %%): %.0f (%.2f%%)", length(igraph::E(g)), igraph::edge_density(g)*100),
      "\n ", sprintf("         N components: %.0f", igraph::components(g)$no),
      "\n",
      sep = ""
    )
  }

  ### summaries
  net_summary(g = phylo_net.ig)      # phylo
  net_summary(g = phylpp_net.ig)     # phylo + PP
  net_summary(g = phylpphh_net.ig)   # phylo + PP + HH


#######################################################################

# (1.1) scale-free analysis ~ compute positive rates among contacts --------
  ### Phylo + PP
  ig = phylpp_net.ig
  ig_names = igraph::V(ig)$name
  ig_names = ig_names[!grepl(pattern = "[[:alpha:]]", ig_names)]
  phylpp_net.ever_pos = spec_data[spec_data$EVENT_ID %in% ig_names, ] |>
    dplyr::group_by(EVENT_ID) |>
    dplyr::mutate(pos_stat = dplyr::case_when(grepl(pattern = "(^negative$)|(^not detected$)", tolower(RESULT_NAME)) ~ 0,
                                              grepl(pattern = "(^positive$)|(^detected$)|(^severe acute respiratory$)", tolower(RESULT_NAME)) ~ 1)) |>
    dplyr::summarise(ever_pos = any(pos_stat == 1)) |>
    dplyr::ungroup() |>
    dplyr::pull(ever_pos)
  phylpp_net.ever_pos = phylpp_net.ever_pos[!is.na(phylpp_net.ever_pos)]


  ### Phylo + PP + HH
  ig = phylpphh_net.ig
  ig_names = igraph::V(ig)$name
  ig_names = ig_names[!grepl(pattern = "[[:alpha:]]", ig_names)]
  phylpphh_net.ever_pos = spec_data[spec_data$EVENT_ID %in% ig_names, ] |>
    dplyr::group_by(EVENT_ID) |>
    dplyr::mutate(pos_stat = dplyr::case_when(grepl(pattern = "(^negative$)|(^not detected$)", tolower(RESULT_NAME)) ~ 0,
                                              grepl(pattern = "(^positive$)|(^detected$)|(^severe acute respiratory$)", tolower(RESULT_NAME)) ~ 1)) |>
    dplyr::summarise(ever_pos = any(pos_stat == 1)) |>
    dplyr::ungroup() |>
    dplyr::pull(ever_pos)
  phylpphh_net.ever_pos = phylpphh_net.ever_pos[!is.na(phylpphh_net.ever_pos)]


  ### print
  cat(
    "Positivity Rate among contacts (non-phylos) \n",
    " * Phylo + PP: ", mean(phylpp_net.ever_pos, na.rm = TRUE), "\n",
    " * Phylo + PP + HH: ", mean(phylpphh_net.ever_pos, na.rm = TRUE), "\n",
    sep = ""
  )

#######################################################################

# (2.0) visuals .............. -----
  ### generic function to visualize networks
  network_vis = function(ig,
                         layout_trans, layout_seed,
                         node_size, phyl.alpha, non_phyl.alpha, edge.alpha,
                         node.color, palette,
                         legend.position, title) {
    set.seed(layout_seed)

    # controls
    node_type1 = "Present in the Phylogenetic-Only network"
    node_type2 = "Included via connections through social networks"
    igraph::V(ig)$Type = ifelse(grepl("phyl", igraph::V(ig)$set), node_type1, node_type2)
    igraph::V(ig)$Type.alpha = ifelse(grepl("phyl", igraph::V(ig)$set), phyl.alpha, non_phyl.alpha)

    # palette
    if (!is.null(palette)) {
      names(palette) = c(node_type1, node_type2)
    }

    # layouts
    layout_coords = igraph::layout_with_fr(graph = ig)
    layout_trans_mat = matrix(layout_trans, byrow = T, nrow = length(igraph::V(ig)), ncol = 2)

    # plot
    vis = GGally::ggnet2(
      # graph & positions
      net = ig,
      mode = layout_coords * layout_trans_mat,
      # edge controls
      edge.alpha = edge.alpha,
      edge.color = "gray65",
      # node controls
      node.alpha = "Type.alpha",
      node.color = node.color,
      node.shape = 19,
      size = node_size,
      # misc.
      palette = palette
    )

    # control aesthetics
    vis = vis +
      ggplot2::theme_bw() +
      ggplot2::labs(title = title, x = "", y = "") +
      ggplot2::theme(title = ggplot2::element_text(face = "bold"),
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     legend.position = legend.position)

    # control legend
    if (legend.position != "none") {
      vis = vis + ggplot2::guides(color = ggplot2::guide_legend(ncol = 1, override.aes = list(size = 3)))
    }

    # output
    vis
  }

  ### visualize level 1 - 3
  vis.L1 = network_vis(
    ig = phylo_net.ig,
    layout_seed = 123, layout_trans = c(1, 1),
    node_size = 1, phyl.alpha = 1, non_phyl.alpha = 0.75, edge.alpha = 0.75,
    node.color = "red3", palette = NULL,
    legend.position = "none", title = ""
  )
  vis.L2 = network_vis(
    ig = phylpp_net.ig,
    layout_seed = 123, layout_trans = c(1, 1),
    node_size = 1, phyl.alpha = 1, non_phyl.alpha = 0.75, edge.alpha = 0.75,
    node.color = "Type", palette = c("red3", "gray50"),
    legend.position = "bottom", title = ""
  )
  vis.L3 = network_vis(
    ig = phylpphh_net.ig,
    layout_seed = 145, layout_trans = c(-1, 1),
    node_size = 1, phyl.alpha = 1, non_phyl.alpha = 0.75, edge.alpha = 0.75,
    node.color = "Type", palette = c("red3", "gray50"),
    legend.position = "bottom", title = ""
  )


  ### all together as 1x3 singal image
  vis.combined.L1to3 = ggpubr::ggarrange(
    vis.L1,
    ggpubr::ggarrange(vis.L2, vis.L3, nrow = 2, ncol = 1, legend = "none"),
    nrow = 1, ncol = 2,
    widths = c(1, 0.8),
    legend = "bottom",
    legend.grob = ggpubr::get_legend(vis.L2)
  )


  # ### save
  # ggplot2::ggsave(
  #   filename = paste(path.save, "/Visual_Phylo_v", today_yyyymmdd, ".png", sep = ""),
  #   plot = vis.L1,
  #   width = 7, height = 7.1,
  #   device = "png",
  #   dpi = "retina"
  # )
  # ggplot2::ggsave(
  #   filename = paste(path.save, "/Visual_PhyloPP_v", today_yyyymmdd, ".png", sep = ""),
  #   plot = vis.L2,
  #   width = 6.8, height = 7.4,
  #   device = "png",
  #   dpi = "retina"
  # )
  # ggplot2::ggsave(
  #   filename = paste(path.save, "/Visual_PhyloPPHH_v", today_yyyymmdd, ".png", sep = ""),
  #   plot = vis.L3,
  #   width = 6.8, height = 7.4,
  #   device = "png",
  #   dpi = "retina"
  # )
  # ggplot2::ggsave(
  #   filename = paste(path.save, "/Visual_threeLayers_combined_v", today_yyyymmdd, ".png", sep = ""),
  #   plot = vis.combined.L1to3,
  #   width = 12, height = 10,
  #   device = "png",
  #   dpi = "retina"
  # )


##################################



# (3.0) scale-free analysis ~ GoF ............. -----
sink(file = paste(path.save, "/ScaleFree_analysis_3layer_network_v", today_yyyymmdd, ".txt", sep = ""))

  ### network data
  pl_net_list = list(phylo_net.ig = phylo_net.ig,
                     phylpp_net.ig = phylpp_net.ig,
                     phylpphh_net.ig = phylpphh_net.ig)

  ### loop for each network
  result_list = list()
  for(i in seq_along(pl_net_list)){
    # current network + observed degree
    net = pl_net_list[[i]]
    my_k = igraph::degree(net) |> as.numeric()

    # fit PL via igraph
    my_fit.PL = igraph::fit_power_law(x = my_k, xmin = NULL)

    # GoF
    my_fit.PL.gof = SFdegree::gof_PL(
      k = my_k,
      kmin = my_fit.PL$xmin,
      alpha = my_fit.PL$alpha,
      Nsim = 2500,
      ss = 123
    )

    # GoF's p-value is proportion of KS_sim > KS_obs
    #   KS_sim = KS statistic from simulated data
    #   KS_obs = KS statistic from actual data
    cat(
      "\n",
      " ~ SUMMARY ~ \n",
      "  * name: ", names(pl_net_list)[i], "\n",
      "  * k min was equal to: ", my_fit.PL$xmin, " (a=", sprintf("%.3f", my_fit.PL$alpha), ") \n",
      "    - thus, ", sum(my_k >= my_fit.PL$xmin), "/", length(my_k), " (", sprintf("%.3f", 100*mean(my_k >= my_fit.PL$xmin)), "%) of the original data were analyzed. \n",
      "  * the GoF p-value was: ", mean(my_fit.PL.gof[, "KS.stat"] > my_fit.PL$KS.stat), "\n",
      " ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \n\n",
      sep = ""
    )

    result_list[[i]] = list(PL_model = my_fit.PL, GoF_stat = my_fit.PL.gof)
  }

sink(file = NULL)
###################################################

# (3.1) scale-free analysis ~ alt models ...... ------
  ### function to conduct the scale-free analysis + compute measures relevant statistics
  test_scale_free = function(net,
                             xmin, force.continuous,
                             exp.inits = matrix(), lnorm.inits = matrix(),
                             PLc.inits = matrix(), weib.inits = matrix(),
                             plot.right_max,
                             jitter = list(jitter = FALSE, jitter_max = NA),
                             x_ticks, y_ticks,
                             y_type = c("pdf", "cdf", "ccdf"),
                             axis_scales = list(xlog = FALSE, ylog = FALSE)) {
    # a list object holding the best fit for each model
    best_fits = list()

    # degree & fit power law
    my_k = as.numeric(igraph::degree(net))
    my_fit.PL = igraph::fit_power_law(x = my_k, xmin = xmin, force.continuous = force.continuous)
    my_fit.PL$par = c(kmin = my_fit.PL$xmin, alpha = my_fit.PL$alpha)
    my_fit.PL$LL = my_fit.PL$logLik
    best_fits$PL = my_fit.PL

    # fit exp & find best fit with corresponding MLE
    my_fit.exp = apply(
      X = exp.inits,
      MARGIN = 1,
      FUN = function(i) {
        SFdegree::fit_exp(k = my_k, kmin = my_fit.PL$xmin, inits = c(log.lambda = i))
      }
    )
    best_fits$exp = my_fit.exp[[which.max(sapply(my_fit.exp, function(x) x$LL))]]

    # fit log-normal & find best fit with corresponding MLE
    my_fit.lnorm = apply(
      X = lnorm.inits,
      MARGIN = 1,
      FUN = function(x) {
        SFdegree::fit_lognorm(k = my_k, kmin = my_fit.PL$xmin, inits = c(mu = x[1], log.sigma = x[2]))
      }
    )
    best_fits$lnorm = my_fit.lnorm[[which.max(sapply(my_fit.lnorm, function(x) x$LL))]]

    # fit power-law with exponential cutoff & find best fit with corresponding MLE
    my_fit.PLcut = apply(
      X = PLc.inits,
      MARGIN = 1,
      FUN = function(x) {
        SFdegree::fit_PLc(k = my_k, kmin = my_fit.PL$xmin, inits = c(log.alpha = x[1], log.lambda = x[2]))
      }
    )
    best_fits$PLcut = my_fit.PLcut[[which.max(sapply(my_fit.PLcut, function(x) x$LL))]]

    # fit weibull & find best fit with corresponding MLE
    my_fit.weib = apply(
      X = weib.inits,
      MARGIN = 1,
      FUN = function(x) {
        SFdegree::fit_weib(k = my_k, kmin = my_fit.PL$xmin, inits = c(log.a = x[1], log.b = x[2]))
      }
    )
    best_fits$weib = my_fit.weib[[which.max(sapply(my_fit.weib, function(x) x$LL))]]


    # MLEs for each model
    my_fit.exp.par = best_fits$exp$par
    my_fit.lnorm.par = best_fits$lnorm$par
    my_fit.PLcut.par = best_fits$PLcut$par
    my_fit.weib.par = best_fits$weib$par

    # compute LLi for each model
    LLi = list(
      PL = SFdegree::PL_LLi(k = my_k,
                            kmin = best_fits$PL$par["kmin"],
                            params = list(alpha = best_fits$PL$par["alpha"]), log = TRUE),
      exp = SFdegree::exp_LLi(k = my_k,
                              kmin = best_fits$PL$par["kmin"],
                              params = list(lambda = best_fits$exp$par["lambda"]), log = TRUE),
      lnorm = SFdegree::lognorm_LLi(k = my_k,
                                    kmin = best_fits$PL$par["kmin"],
                                    params = list(mu = best_fits$lnorm$par["mu"], sigma = best_fits$lnorm$par["sigma"]), log = TRUE),
      PLcut = SFdegree::PLc_LLi(k = my_k,
                                kmin = best_fits$PL$par["kmin"],
                                params = list(alpha = best_fits$PLcut$par["alpha"], lambda = best_fits$PLcut$par["lambda"]), log = TRUE),
      weib = SFdegree::weib_LLi(k = my_k,
                                kmin = best_fits$PL$par["kmin"],
                                params = list(a = best_fits$weib$par["a"], b = best_fits$weib$par["b"]), log = TRUE)
    )

    # 5 choose 2 LRTs
    LRT_mat.R = matrix(NA, nrow = 5, ncol = 5)
    rownames(LRT_mat.R) = colnames(LRT_mat.R) = names(best_fits)
    diag(LRT_mat.R) = 0

    LRT_mat.P = matrix(NA, nrow = 5, ncol = 5)
    rownames(LRT_mat.P) = colnames(LRT_mat.P) = names(best_fits)
    diag(LRT_mat.P) = "-"
    for (i in rownames(LRT_mat.R)) {
      for (j in colnames(LRT_mat.R)) {
        if (i != j) {
          lrt = SFdegree::LRTest(LLi1 = LLi[[i]], LLi2 = LLi[[j]])

          LRT_mat.R[i, j] = lrt$Rstat
          LRT_mat.P[i, j] = ifelse(lrt$pval < 0.001, "<0.001", sprintf("%.4f", lrt$pval))
        }
      }
    }

    # the actual degree distribution data
    set.seed(123)
    my_k.right = best_fits$exp$data # vector of the right-tail data, not the original
    k_range = seq(range(my_k.right)[1], range(my_k.right)[2], 1)
    k_freq = as.numeric(table(factor(my_k.right, levels = k_range)) / length(my_k.right))

    # a data frame containing the densitiy for each model
    density_data = data.frame(x = k_range,
                              our_data = k_freq,
                              PL = SFdegree::PL_LLi(k = k_range,
                                                    kmin = best_fits$PL$par["kmin"],
                                                    params = list(alpha = as.numeric(best_fits$PL$alpha)),
                                                    log = FALSE),
                              Exp = SFdegree::exp_LLi(k = k_range,
                                                      kmin = best_fits$PL$par["kmin"],
                                                      params = list(lambda = best_fits$exp$par["lambda"]),
                                                      log = FALSE),
                              LogNorm = SFdegree::lognorm_LLi(k = k_range,
                                                              kmin = best_fits$PL$xmin,
                                                              params = list(mu = best_fits$lnorm$par["mu"],
                                                                            sigma = best_fits$lnorm$par["sigma"]),
                                                              log = FALSE),
                              PLCut = SFdegree::PLc_LLi(k = k_range,
                                                        kmin = best_fits$PL$par["kmin"],
                                                        params = list(alpha = best_fits$PLcut$par["alpha"],
                                                                      lambda = best_fits$PLcut$par["lambda"]),
                                                        log = FALSE),
                              Weib = SFdegree::weib_LLi(k = k_range,
                                                        kmin = best_fits$PL$par["kmin"],
                                                        params = list(a = best_fits$weib$par["a"],
                                                                      b = best_fits$weib$par["b"]),
                                                        log = FALSE))

    # formats density data frame
    density_data = density_data |>
      tidyr::pivot_longer(cols = -c(x), names_to = "Models", values_to = "density")
    if (jitter$jitter) {
      density_data = density_data |>
        dplyr::group_by(Models) |>
        dplyr::mutate(
          jit = runif(n = 1, min = 0, max = jitter$jitter_max),
          density = density + jit
        ) |>
        dplyr::ungroup()
    }
    density_data = density_data |>
      dplyr::group_by(Models) |>
      dplyr::mutate(
        cumulative_prob = cumsum(density),
        cumulative_prob = cumulative_prob / max(cumulative_prob)
      ) |>
      dplyr::ungroup()


    # data frame to be plotted
    plot_data = density_data |> dplyr::filter(x <= plot.right_max)
    if (y_type == "pdf") {
      plot_data$y = plot_data$density
    }
    if (y_type == "cdf") {
      plot_data$y = plot_data$cumulative_prob
    }
    if (y_type == "ccdf") {
      plot_data$y = 1 - plot_data$cumulative_prob
    }

    # plot models
    p = ggplot2::ggplot() +
      ggplot2::geom_path(data = plot_data |> dplyr::filter(Models != "our_data"),
      # ggplot2::geom_path(data = plot_data |> dplyr::filter(!Models %in% c("our_data", "Weib")),
                         ggplot2::aes(x = x, y = y, group = Models, color = Models),
                         linewidth = 1, linetype = "solid", alpha = 0.7)
    # plot empirical data
    p = p +
      ggplot2::geom_point(data = plot_data |> dplyr::filter(Models == "our_data") ,
                          ggplot2::aes(x = x, y = y),
                          size = 1.5, alpha = 0.9, color = 'gray30') +
      ggplot2::geom_path(data = plot_data |> dplyr::filter(Models == "our_data") ,
                         ggplot2::aes(x = x, y = y),
                         linewidth = 0.7, linetype = 'dashed', alpha = 1, color = 'black')

    # axis scaling controls
    if (axis_scales$xlog) {
      p = p + ggplot2::scale_x_log10(breaks = x_ticks,
                                     limits = range(x_ticks),
                                     labels = scales::trans_format("log10", scales::math_format(10^.x)))
    } else {
      p = p + ggplot2::scale_x_continuous(breaks = x_ticks,
                                          limits = range(x_ticks))
    }
    if (axis_scales$ylog) {
      p = p + ggplot2::scale_y_log10(breaks = y_ticks,
                                     limits = range(y_ticks),
                                     labels = scales::trans_format("log10", scales::math_format(10^.x)))
    }

    # other plot controls
    p = p +
      ggplot2::labs(title = "", x = "", y = "") +
      ggplot2::theme_bw() +
      ggplot2::theme(title = ggplot2::element_text(face = "bold"),
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     legend.position = "bottom")


    # output
    out = list(
      best_fits = best_fits,
      n_par = lapply(best_fits, FUN = function(x) length(x$par)),
      MLE = lapply(best_fits, FUN = function(x) x$par),
      LL = lapply(best_fits, FUN = function(x) x$LL),
      AIC = lapply(best_fits, FUN = function(x) 2*length(x$par) - 2*x$LL),
      LRT = list(pval = LRT_mat.P, Rstat = LRT_mat.R),
      plot = p,
      data.right_tail = my_k.right
    )
    out
  }

  ### test Phylo Only
  phylo_net.testSF = test_scale_free(
    # network
    net = phylo_net.ig,
    # Power-law model controls
    xmin = NULL, force.continuous = FALSE,
    # initial values for alternative models
    exp.inits = matrix(c(-3.5, -0.8, 0.8, 3.5), byrow = T, nrow = 4),
    lnorm.inits = as.matrix(expand.grid(c(-10, -1, 1, 10), c(-6, -1, 1, 6)) |> unname()),
    PLc.inits = as.matrix(expand.grid(c(-6, -2, 2, 6), c(-4, -1, 1, 4)) |> unname()),
    weib.inits = as.matrix(expand.grid(c(-10, -1, 0.5, 2), c(-7, -1, 2, 6)) |> unname()),
    # plot controls
    plot.right_max = 9999,
    jitter = list(jitter = F, jitter_max = 99),
    x_ticks = 10^seq(0, 0.8, 0.2),
    y_ticks = 10^seq(-2.5, 0, 0.5),
    y_type = "ccdf",
    axis_scales = list(xlog = T, ylog = T)
  )
  
  ### test Phylo + PP
  phylpp_net.testSF = test_scale_free(
    # network
    net = phylpp_net.ig,
    # Power-law model controls
    xmin = NULL, force.continuous = FALSE,
    # initial values for alternative models
    exp.inits = matrix(c(-3.5, -0.8, 0.8, 3.5), byrow = T, nrow = 4),
    lnorm.inits = as.matrix(expand.grid(c(-10, -1, 1, 10), c(-6, -1, 1, 6)) |> unname()),
    PLc.inits = as.matrix(expand.grid(c(-6, -2, 2, 6), c(-4, -1, 1, 4)) |> unname()),
    weib.inits = as.matrix(expand.grid(c(-10, -1, 0.5, 2), c(-7, -1, 2, 6)) |> unname()),
    # plot controls
    plot.right_max = 9999,
    jitter = list(jitter = F, jitter_max = 0.0005),
    x_ticks = 10^seq(0.4, 1.9, 0.3),
    y_ticks = 10^seq(-5, 0, 1),
    y_type = "ccdf",
    axis_scales = list(xlog = T, ylog = T)
  )


  ### test Phylo + PP + HH
  phylpphh_net.testSF = test_scale_free(
    # network
    net = phylpphh_net.ig,
    # Power-law model controls
    xmin = NULL, force.continuous = FALSE,
    # initial values for alternative models
    exp.inits = matrix(c(-3.5, -0.8, 0.8, 3.5), byrow = T, nrow = 4),
    lnorm.inits = as.matrix(expand.grid(c(-10, -1, 1, 10), c(-6, -1, 1, 6)) |> unname()),
    PLc.inits = as.matrix(expand.grid(c(-6, -2, 2, 6), c(-4, -1, 1, 4)) |> unname()),
    weib.inits = as.matrix(expand.grid(c(-10, -1, 0.5, 2), c(-7, -1, 2, 6)) |> unname()),
    # plot controls
    plot.right_max = 9999,
    jitter = list(jitter = F, jitter_max = 0.002),
    x_ticks = 10^seq(0.3, 1.9, 0.4),
    y_ticks = 10^seq(-5, 0, 1),
    y_type = "ccdf",
    axis_scales = list(xlog = T, ylog = T)
  )



###################################################

# (3.2) scale-free analysis ~ plot densities .. --------------
### unrestricted Phylo + PP + HH plot
phylpphh_net.testSF.v2 = test_scale_free(
  # network
  net = phylpphh_net.ig,
  # Power-law model controls
  xmin = NULL, force.continuous = FALSE,
  # initial values for alternative models
  exp.inits = matrix(c(-3.5, -0.8, 0.8, 3.5), byrow = T, nrow = 4),
  lnorm.inits = as.matrix(expand.grid(c(-10, -1, 1, 10), c(-6, -1, 1, 6)) |> unname()),
  PLc.inits = as.matrix(expand.grid(c(-6, -2, 2, 6), c(-4, -1, 1, 4)) |> unname()),
  weib.inits = as.matrix(expand.grid(c(-10, -1, 0.5, 2), c(-7, -1, 2, 6)) |> unname()),
  # plot controls
  plot.right_max = 9999,
  jitter = list(jitter = F, jitter_max = 0.002),
  x_ticks = 10^seq(0.3, 1.9, 0.4),
  y_ticks = 10^seq(-5, 0, 1),
  y_type = "ccdf",
  axis_scales = list(xlog = T, ylog = T)
)


### altogether as 2x2 single figure
final.density_plot.all3 = ggpubr::ggarrange(
  phylo_net.testSF$plot + ggplot2::labs(title = "Phylo Only", x = "", y = ""),
  ggpubr::ggarrange(phylpp_net.testSF$plot + ggplot2::labs(title = "Phylo + PP", x = "", y = ""),
                    phylpphh_net.testSF$plot + ggplot2::labs(title = "Phylo + PP + HH", x = "", y = ""),
                    ncol = 1, nrow = 2,
                    legend = "none"),
  nrow = 1, ncol = 2,
  widths = c(1, 0.8),
  legend = "top",
  legend.grob = ggpubr::get_legend(phylpphh_net.testSF$plot)
)
final.density_plot.all3 = final.density_plot.all3 |> ggpubr::annotate_figure(left = "Probability", bottom = "degree")


### compare restricted X range vs. non restricted
compare_unrestricted_x = ggpubr::ggarrange(
  phylpphh_net.testSF$plot + ggplot2::labs(title = "Phylo + PP + HH", x = "", y = ""),
  phylpphh_net.testSF.v2$plot + ggplot2::labs(title = "Phylo + PP + HH", x = "", y = ""),
  common.legend = TRUE,
  legend = "top",
  nrow = 1, ncol = 2
)
compare_unrestricted_x = compare_unrestricted_x |> ggpubr::annotate_figure(left = "Probability", bottom = "degree", )


# ### save
# # combined
# ggplot2::ggsave(
#   plot = final.density_plot.all3,
#   filename = paste(path.save, "/Combined_DegreeDensity_all3_Plot.png", sep = ""),
#   width = 11, height = 9,
#   device = "png",
#   dpi = "retina"
# )
# # individuals
# w = 7.7; h = 7
# ggsave(
#   plot = phylo.density_plot + labs(title = "Phylo Only", x = "", y = "") + theme(title = element_text(face = "bold")),
#   filename = paste(path.save, "/Phylo_DegreeDensity_Plot.png", sep = ""),
#   device = "png",
#   dpi = "retina",
#   width = w, height = h
# )
# ggsave(
#   plot = phylpp.density_plot + labs(title = "Phylo + PP", x = "", y = "") + theme(title = element_text(face = "bold")),
#   filename = paste(path.save, "/PhylPP_DegreeDensity_Plot.png", sep = ""),
#   device = "png",
#   dpi = "retina",
#   width = w, height = h
# )
# ggsave(
#   plot = phylpphh.density_plot + labs(title = "Phylo + PP + HH", x = "", y = "") + theme(title = element_text(face = "bold")),
#   filename = paste(path.save, "/PhylPPHH_DegreeDensity_Plot.png", sep = ""),
#   device = "png",
#   dpi = "retina",
#   width = w, height = h
# )



###################################################

# (3.3) scale-free analysis ~ supplements (GoF) ------
### Phylo + PP network
net = phylpp_net.ig
my_k = as.numeric(igraph::degree(net))
my_k = my_k[my_k >= 6]

my_fit.PL = igraph::fit_power_law(x = my_k, xmin = 6)

# GoF
my_fit.PL.gof = SFdegree::gof_PL(
  k = my_k,
  kmin = my_fit.PL$xmin,
  alpha = my_fit.PL$alpha,
  Nsim = 250,
  ss = 123
)

cat(
  "\n",
  " ~ SUMMARY ~ \n",
  "  * name: ", "Phylo + PP", "\n",
  "  * k min was equal to: ", my_fit.PL$xmin, " (a=", sprintf("%.3f", my_fit.PL$alpha), ") \n",
  "    - thus, ", sum(my_k >= my_fit.PL$xmin), "/", length(my_k), " (", sprintf("%.3f", 100*mean(my_k >= my_fit.PL$xmin)), "%) of the original data were analyzed. \n",
  "  * the GoF p-value was: ", mean(my_fit.PL.gof[, "KS.stat"] > my_fit.PL$KS.stat), "\n",
  " ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \n\n",
  sep = ""
)


### Phylo + PP + HH network
net = phylpphh_net.ig
my_k = as.numeric(igraph::degree(net))
my_k = my_k[my_k >= 7]

my_fit.PL = igraph::fit_power_law(x = my_k, xmin = 7)

# GoF
my_fit.PL.gof = SFdegree::gof_PL(
  k = my_k,
  kmin = my_fit.PL$xmin,
  alpha = my_fit.PL$alpha,
  Nsim = 250,
  ss = 123
)

cat(
  "\n",
  " ~ SUMMARY ~ \n",
  "  * name: ", "Phylo + PP", "\n",
  "  * k min was equal to: ", my_fit.PL$xmin, " (a=", sprintf("%.3f", my_fit.PL$alpha), ") \n",
  "    - thus, ", sum(my_k >= my_fit.PL$xmin), "/", length(my_k), " (", sprintf("%.3f", 100*mean(my_k >= my_fit.PL$xmin)), "%) of the original data were analyzed. \n",
  "  * the GoF p-value was: ", mean(my_fit.PL.gof[, "KS.stat"] > my_fit.PL$KS.stat), "\n",
  " ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \n\n",
  sep = ""
)

###################################################

# (3.4) scale-free analysis ~ combined figure . -----
"Pre-req to run this section is that section (2), (3.1), (3.2) must be ran beforehand"

### row 1
plot3x2.r1 = ggpubr::ggarrange(
  vis.L1, phylo_net.testSF$plot,
  legend = 'none',
  ncol = 2, nrow = 1
)

### row 2
plot3x2.r2 = ggpubr::ggarrange(
  vis.L2, phylpp_net.testSF$plot,
  legend = 'none',
  ncol = 2, nrow = 1
)

### row 3
plot3x2.r3 = ggpubr::ggarrange(
  vis.L3, phylpphh_net.testSF$plot,
  legend = 'none',
  ncol = 2, nrow = 1
)

### all together (3x2)
plot3x2 = ggpubr::ggarrange(
  plot3x2.r1, plot3x2.r2, plot3x2.r3,
  labels = list("Phylo Only", "Phylo + PP", "Phylo + PP + HH"),
  hjust = c(-4.4, -4.4, -2.8),
  vjust = c(1.2, 1.2, 1.2),
  ncol = 1, nrow = 3
)

# ### save
# ggplot2::ggsave(
#   filename = paste(path.save, "/COVID-19 network figures 3x2 v", today_yyyymmdd, ".png", sep = ""),
#   plot = plot3x2,
#   device = "png", dpi = "retina",
#   height = 13, width = 9.75,
# )

###################################################


