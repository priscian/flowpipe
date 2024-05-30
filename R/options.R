### For later expansion

.flowpipe <- structure(
  list(
    params = NULL
  )) %>% keystone::add_class("package_var")


# `$<-.package_var` <- function(x, name, value)
# {
#   # print("`$<-.package_var`")
#   current <- x
#   current[[name]] <- value
#   assign(".flowpipe", current, envir = asNamespace("flowpipe"))
# }


#' @export
flowpipe_params <- keystone::package_params %>%
  `formals<-`(value = formals(.) %>%
    `$<-`("__NAMESPACE__", "flowpipe") %>%
    `$<-`("__VARNAME__", ".flowpipe")) %>%
  `environment<-`(asNamespace("flowpipe"))
