Person <- R6::R6Class("Person",
                      public = list(
#' @description 
#' Constructor of the 'Person' class
#' @param name a string that defines this person's name
#' @param id an integer that defines this person's id
#' @return A new 'Person' object
                        initialize = function(val){private$value <- val},
                        get_value = function(){return(private$name)},
                        
                        set_value = function(){return(private$id)}
                      )
)
  