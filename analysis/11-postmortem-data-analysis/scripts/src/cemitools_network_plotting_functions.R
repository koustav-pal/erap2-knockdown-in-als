setGeneric('plot_interactions', function(cem, ...) {
    standardGeneric('plot_interactions')
})

#' @rdname plot_interactions
setMethod('plot_interactions', signature('CEMiTool'),
    function(cem, n=10, label_genes, ...) {
        if(length(unique(cem@module$modules)) == 0){
            stop("No modules in CEMiTool object! Did you run find_modules()?")
        }
        if(length(interactions_data(cem)) == 0){
            stop("No interactions information! Did you run interactions_data()?")
        }
        #cem <- get_args(cem, vars=mget(ls()))
        mod_cols <- mod_colors(cem)
        mod_names <- names(cem@interactions)
        mod_names <- mod_names[which(mod_names!="Not.Correlated")]
        hubs <- lapply(get_hubs(cem), names)
        zero_ints <- character()
        zero_ints <- lapply(names(cem@interactions), function(mod){
                        degree <- igraph::degree(cem@interactions[[mod]], normalized=FALSE)
                        if(length(degree) == 0) {
                            zero_ints <- append(zero_ints, mod)
                        }
                     })
        zero_ints <- unlist(zero_ints)
        if(!is.null(zero_ints)){
            mod_names <- mod_names[which(!(mod_names %in% zero_ints))]
        }
        if(length(mod_names) == 0){
            warning("There are no interactions in the given modules. Please check interactions file.")
            return(cem)
        }
        res <- lapply(mod_names, function(name){
                   plot_interaction(ig_obj=cem@interactions[[name]],
                                    n=n, color=mod_cols[name], name=name,
                                    mod_genes=module_genes(cem, name)$genes,
                                    label_genes = label_genes, 
                                    coexp_hubs=hubs[[name]])
               })
        names(res) <- mod_names
        mod_names_ordered <- mod_names[order(as.numeric(stringr::str_extract(mod_names, "\\d+")))]
        cem@interaction_plot <- res[mod_names_ordered]
        return(cem)
    })

plot_interaction <- function(ig_obj, n, color, name, mod_genes, label_genes = NA, coexp_hubs){
    degrees <- igraph::degree(ig_obj, normalized=FALSE)
    ig_obj <- igraph::set_vertex_attr(ig_obj, "degree", value = degrees)
    max_n <- min(n, length(degrees))
    net_obj <- intergraph::asNetwork(ig_obj)
    m <- network::as.matrix.network.adjacency(net_obj) # get sociomatrix
    # get coordinates from Fruchterman and Reingold's force-directed placement algorithm.
    plotcord <- data.frame(sna::gplot.layout.fruchtermanreingold(m, NULL))
    # or get it them from Kamada-Kawai's algorithm:
    # plotcord <- data.frame(sna::gplot.layout.kamadakawai(m, NULL))
    colnames(plotcord) <- c("X1","X2")
    edglist <- network::as.matrix.network.edgelist(net_obj)
    edges <- data.frame(plotcord[edglist[,1],], plotcord[edglist[,2],])
    plotcord$vertex.names <- as.factor(network::get.vertex.attribute(net_obj, "vertex.names"))
    plotcord$Degree <- network::get.vertex.attribute(net_obj, "degree")
    plotcord[, "shouldLabel"] <- FALSE
    plotcord[, "Hub"] <- ""
    int_hubs <- names(sort(degrees, decreasing=TRUE))[1:max_n]
    int_bool <- plotcord[, "vertex.names"] %in% int_hubs
    plotcord[which(int_bool), "Hub"] <- "Interaction"
    sel_vertex <- int_hubs
    if(!missing(coexp_hubs)){
        coexp_bool <- plotcord[, "vertex.names"] %in% coexp_hubs
        coexp_and_int <- coexp_bool & int_bool
        plotcord[which(coexp_bool), "Hub"] <- "Co-expression"
        plotcord[which(coexp_and_int), "Hub"] <- "Co-expression + Interaction"
        sel_vertex <- c(sel_vertex, coexp_hubs)
    }

    colnames(edges) <-  c("X1","Y1","X2","Y2")
    #edges$midX  <- (edges$X1 + edges$X2) / 2
    #edges$midY  <- (edges$Y1 + edges$Y2) / 2
    plotcord[which(plotcord[, "vertex.names"] %in% sel_vertex), "shouldLabel"] <- TRUE
    if(all(!is.na(label_genes))){
        plotcord[which(plotcord[, "vertex.names"] %in% label_genes), "shouldLabel"] <- TRUE    
    }
    plotcord$Degree_cut <- cut(plotcord$Degree, breaks=3, labels=FALSE)
    plotcord$in_mod <- TRUE
    #mod_genes <- cem@module[cem@module$modules==name,]$genes
    not_in <- setdiff(plotcord[,"vertex.names"], mod_genes)
    plotcord[which(plotcord[, "vertex.names"] %in% not_in), "in_mod"] <- FALSE

    pl <- ggplot(plotcord)  +
        geom_segment(data=edges, aes_(x=~X1, y=~Y1, xend=~X2, yend=~Y2),
                     size = 0.5, alpha=0.5, colour="#DDDDDD") +
        geom_point(aes_(x=~X1, y=~X2, size=~Degree, alpha=~Degree), color=color) +
        geom_label_repel(aes_(x=~X1, y=~X2, label=~vertex.names, color=~Hub),
                         box.padding=unit(1, "lines"),
                         data=function(x){x[x$shouldLabel, ]}) +
        scale_colour_manual(values=c("Co-expression" = "#005E87",
                                     "Interaction" = "#540814",
                                     "Co-expression + Interaction" = "#736E0B")) +
        labs(title=name) +
        ggplot2::theme_bw(base_size = 12, base_family = "") +
        ggplot2::theme(axis.text = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.title = ggplot2::element_blank(),
                       legend.key = ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill = "white",
                                                                colour = NA),
                       panel.border = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank())

    return(pl)
}
