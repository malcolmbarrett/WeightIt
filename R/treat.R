#Treatment class

#The class every processed treatment carries. It is cobalt's (see
#`cobalt::treat-class`), and so is the `[` method that preserves the attributes across
#subsetting: cobalt registers it and WeightIt does not, so there is one method for one
#class and no chance of the two packages overwriting each other's.
#
#`cobalt.treat` is the class the method is currently registered on, and `treat` is the
#shared contract. Both are set, so the object is indistinguishable from one cobalt
#processed itself and dispatch finds the method whichever of the two names cobalt
#registers it under. It goes first because a multi-category treatment is a factor
#underneath, and `[.factor` would otherwise win and drop every attribute.
.treat_classes <- c("cobalt.treat", "treat")

.set_treat_class <- function(x) {
  class(x) <- unique(c(.treat_classes, class(x)))

  x
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

  .set_treat_class(x)
}
