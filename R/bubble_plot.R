bubble_plot <- function(object, mod, group = NULL, xlab = "Moderator", ylab = "Effect size", N = "none",
                        alpha = 0.5, cb = TRUE, k = TRUE, g = FALSE,
                        est.lwd = 1, ci.lwd = 0.5, pi.lwd = 0.5,
                        est.col = "black", ci.col = "black", pi.col = "black",
                        legend.pos = c("top.left", "top.right",
                                       "bottom.right", "bottom.left",
                                       "top.out", "bottom.out",
                                       "none"),
                        k.pos = c("top.right", "top.left",
                                  "bottom.right", "bottom.left",
                                  "none"),
                        condition.nrow = 2,
                         #condition.lab = "Condition",
                        weights = "prop", by = NULL, at = NULL)
  {
  legend.pos <- match.arg(NULL, choices = legend.pos)
  k.pos <- match.arg(NULL, choices = k.pos)
  #facet <- match.arg(NULL, choices = facet)

  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?bubble_plot")
  }

  if(is.numeric(by)){
   k = FALSE
   g = FALSE
  }


  if(any(class(object) %in% c("robust.rma", "rma.mv", "rma", "rma.uni"))){

    if(mod != "1"){
      results <-  orchaRd::mod_results(object, mod, group,
                                       by = by, at = at, weights = weights)
    } else {
      results <-  orchaRd::mod_results(object, mod = "1", group,
                                       by = by, at = at, weights = weights)
    }
  }

  if(any(class(object) %in% c("orchard"))) {
    results <- object
  }

  mod_table <- results$mod_table

  data_trim <- results$data

  data_trim$scale <- (1/sqrt(data_trim[,"vi"]))
  legend <- "Precision (1/SE)"

  if(any(N != "none")){
    data_trim$scale <- data_trim$N
    legend <- paste0("Sample Size ($\\textit{N}$)") # we want to use italic
  }

  label <- xlab
  # if(transfm == "tanh"){
  #   cols <- sapply(mod_table, is.numeric)
  #   mod_table[,cols] <- Zr_to_r(mod_table[,cols])
  #   data_trim$yi <- Zr_to_r(data_trim$yi)
  #   label <- xlab
  # }else{
  #   label <- xlab
  # }

  if(is.null(data_trim$condition) == TRUE){

  # the number of effect sizes
  effect_num <- nrow(data_trim)

  # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
  group_num <- length(unique(data_trim$stdy))

  dat_text <- data.frame(K = effect_num, G = group_num)

  }else{

  # making sure factor names match
  data_trim$condition <- factor(data_trim$condition, levels = mod_table$condition, labels = mod_table$condition)

  effect_num <- as.vector(by(data_trim, data_trim[,"condition"], function(x) base::length(x[,"yi"])))

  # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
  #group_num <- c(2,4)
  group_num <- as.vector(by(data_trim, data_trim[,"condition"], function(x) base::length(base::unique(x[,"stdy"]))))

  dat_text <- data.frame(K = effect_num, G = group_num, condition = as.vector(base::levels(data_trim$condition)))
  }
  # the number of groups in a moderator & data points
  #group_no <- length(unique(mod_table[, "name"]))

  #data_no <- nrow(data)

  # # colour blind friendly colours with grey
  # cbpl <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

  if(is.null(data_trim$condition) == TRUE){
   plot <-ggplot2::ggplot() +
    # putting bubbles
     ggplot2::geom_point(data = data_trim, ggplot2::aes(x = moderator, y = yi, size = scale), shape = 21, alpha = alpha, fill = "grey90" ) +
    # prediction interval
     ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerPR), method =  "loess", formula = y~x, se = FALSE, lty =  "dotted", lwd = pi.lwd, colour = pi.col) +
     ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperPR), method =  "loess", formula = y~x, se = FALSE, lty = "dotted", lwd = pi.lwd, colour = pi.col) +
     # confidence interval
     ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerCL), method =  "loess", formula = y~x, se = FALSE,lty = "dashed", lwd = ci.lwd, colour = ci.col) +
     ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperCL), method =  "loess", formula = y~x, se = FALSE, lty ="dashed", lwd = ci.lwd, colour = ci.col) +
     # main line
     ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = estimate), method =  "loess", formula = y~x, se = FALSE, lwd = est.lwd, colour = est.col) +
    #facet_grid(rows = vars(condition)) +
     ggplot2::labs(x = xlab, y = ylab, size = legend, parse = TRUE) +
     ggplot2::guides(fill = "none", colour = "none") +
    # themes
     ggplot2::theme_bw() +
    #theme(legend.position= c(1, 1), legend.justification = c(1, 1)) +
     ggplot2::theme(legend.direction="horizontal") +
    #theme(legend.background = element_rect(fill = "white", colour = "black")) +
     ggplot2::theme(legend.background = ggplot2::element_blank()) +
     ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black", hjust = 0.5, angle = 90))
   } else if(is.character(data_trim$condition) == TRUE || is.factor(data_trim$condition) == TRUE){

      plot <-ggplot2::ggplot() +
      # putting bubbles
      ggplot2::geom_point(data = data_trim, ggplot2::aes(x = moderator, y = yi, size = scale, fill = condition), shape = 21, alpha = alpha) +
      # prediction interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerPR), method =  "loess", formula = y~x, se = FALSE, lty =  "dotted", lwd = pi.lwd, colour = pi.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperPR), method =  "loess", formula = y~x,se = FALSE, lty = "dotted", lwd = pi.lwd, colour = pi.col) +
      # confidence interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerCL), method =  "loess", formula = y~x,se = FALSE,lty = "dashed", lwd = ci.lwd, colour = ci.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperCL), method =  "loess", formula = y~x,se = FALSE, lty ="dashed", lwd = ci.lwd, colour = ci.col) +
      # main line
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = estimate), method =  "loess", formula = y~x, se = FALSE, lwd = est.lwd, colour = est.col) +
      ggplot2::facet_wrap(ggplot2::vars(condition), nrow = condition.nrow) +
      ggplot2::labs(x = xlab, y = ylab, size = legend, parse = TRUE) +
      ggplot2::guides(fill = "none", colour = "none") +
      # themses
      ggplot2::theme_bw() +
      #theme(legend.position= c(1, 1), legend.justification = c(1, 1)) +
      ggplot2::theme(legend.direction="horizontal") +
      #theme(legend.background = element_rect(fill = "white", colour = "black")) +
      ggplot2::theme(legend.background = ggplot2::element_blank()) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black", hjust = 0.5, angle = 90))

    # if(facet == "rows"){
    #   plot <- plot + facet_grid(rows = vars(condition))
    # } else{
    #   plot <- plot + facet_grid(cols = vars(condition))
    # }


  } else{
    plot <-ggplot2::ggplot() +
      # putting bubbles
      #geom_point(data = data, aes(x = moderator, y = yi, size = scale), shape = 21, alpha = alpha, fill = "grey90" ) +
      # prediction interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerPR), method =  "loess", formula = y~x, se = FALSE, lty =  "dotted", lwd = pi.lwd, colour = pi.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperPR), method =  "loess", formula = y~x,se = FALSE, lty = "dotted", lwd = pi.lwd, colour = pi.col) +
      # confidence interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerCL), method =  "loess", formula = y~x,se = FALSE,lty = "dashed", lwd = ci.lwd, colour = ci.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperCL), method =  "loess", formula = y~x,se = FALSE, lty ="dashed", lwd = ci.lwd, colour = ci.col) +
      # main line
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = estimate), method =  "loess", formula = y~x, se = FALSE, lwd = est.lwd, colour = est.col) +
      ggplot2::facet_wrap(ggplot2::vars(condition), nrow = condition.nrow) +
      ggplot2::labs(x = xlab, y = ylab, size = legend, parse = TRUE) +
      ggplot2::guides(fill = "none", colour = "none") +
      # themses
      ggplot2::theme_bw() # +
      #theme(legend.position= c(1, 1), legend.justification = c(1, 1)) +
      # theme(legend.direction="horizontal") +
      # #theme(legend.background = element_rect(fill = "white", colour = "black")) +
      # theme(legend.background = element_blank()) +
      # theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, angle = 90))
  }

  # adding legend
  if(legend.pos == "bottom.right"){
    plot <- plot + ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
  } else if ( legend.pos == "bottom.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 0), legend.justification = c(0, 0))
  } else if ( legend.pos == "top.right") {
    plot <- plot + ggplot2::theme(legend.position= c(1, 1), legend.justification = c(1, 1))
  } else if (legend.pos == "top.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1))
  } else if (legend.pos == "top.out") {
    plot <- plot + ggplot2::theme(legend.position="top")
  } else if (legend.pos == "bottom.out") {
    plot <- plot + ggplot2::theme(legend.position="bottom")
  } else if (legend.pos == "none") {
    plot <- plot + ggplot2::theme(legend.position="none")
  }

  # putting k and g in
  # c("top.right", "top.left", "bottom.right", "bottom.left","none")
  if(k == TRUE && g == FALSE && k.pos == "top.right"){
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                        mapping = ggplot2::aes(x = Inf, y = Inf),
                        label =  paste("italic(k)==", dat_text$K),
                        parse = TRUE,
                        hjust   = 2,
                        vjust   = 2.5
                        )

  } else if(k == TRUE && g == FALSE && k.pos == "top.left") {
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = Inf),
                         label =  paste("italic(k)==", dat_text$K),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = 2.5
      )
  } else if(k == TRUE && g == FALSE && k.pos == "bottom.right") {
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = Inf, y = -Inf),
                         label =  paste("italic(k)==", rev(dat_text$K)),
                         parse = TRUE,
                         hjust   = 2,
                         vjust   = -1.5
      )
  } else if (k == TRUE && g == FALSE && k.pos == "bottom.left"){
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = -Inf),
                         label =  paste("italic(k)==", dat_text$K),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = -1.5
      )
    # below get g ----

  } else if (k == TRUE && g == TRUE && k.pos == "top.right"){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                                   mapping = ggplot2::aes(x = Inf, y = Inf),
                                   label =  paste("italic(k)==", dat_text$K,
                                                         "~","(", dat_text$G, ")"),
                                   parse = TRUE,
                                   hjust   = 1.5,
                                   vjust   = 2)

  } else if (k == TRUE && g == TRUE && k.pos == "top.left"){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = Inf),
                         label =  paste("italic(k)==", dat_text$K,
                                        "~","(", dat_text$G, ")"),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = 2)
  } else if (k == TRUE && g == TRUE && k.pos == "bottom.right"){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = Inf, y = -Inf),
                         label =  paste("italic(k)==", dat_text$K,
                                        "~","(", dat_text$G, ")"),
                         parse = TRUE,
                         hjust   = 1.5,
                         vjust   = -0.5)
  } else if (k == TRUE && g == TRUE && k.pos == "bottom.left"){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = -Inf),
                         label =  paste("italic(k)==", dat_text$K,
                                        "~","(", dat_text$G, ")"),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = -0.5)
  }

  # # putting colors in
  # if(cb == TRUE){
  #   plot <- plot +
  #     ggplot2::scale_fill_manual(values=cbpl) +
  #     ggplot2::scale_colour_manual(values=cbpl)
  # }

  return(plot)
}
