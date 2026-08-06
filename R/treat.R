#Treatment class
#' @exportS3Method `[` treat
`[.treat` <- function(x, ..., value) {
  y <- NextMethod("[")
  attr(y, "treat.type") <- .attr(x, "treat.type")
  attr(y, "treat.name") <- .attr(x, "treat.name")
  attr(y, "treated") <- .attr(x, "treated")
  attr(y, "control") <- .attr(x, "control")

  class(y) <- class(x)

  y
}

as.treat <- function(x, process = NULL, censoring = NULL) {
  if (is_null(process)) {
    process <- !inherits(x, "treat")
  }

  arg::arg_flag(process)

  if (process || !has_treat_type(x)) {
    #Multi-category treatments are passed through `factor()`, here and inside
    #`assign_treat_type()`, which drops every attribute; the treatment's name is
    #carried across by hand so that it survives processing as it does for the other
    #treatment types. `weightitMSM()` names its `treat.list` from it.
    treat.name <- .attr(x, "treat.name")

    x <- assign_treat_type(x, censoring = censoring)
    treat.type <- get_treat_type(x)

    if (treat.type %in% c("multinomial", "multi-category")) {
      x <- assign_treat_type(factor(x))
    }

    attr(x, "treat.name") <- treat.name
  }

  class(x) <- unique(c("treat", class(x)))

  x
}
