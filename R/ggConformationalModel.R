










#' Title
#'
#' @param model
#' @param id_states
#' @param delta_z
#' @param dz_x
#' @param dz_y
#'
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
g_ConformationalStates <- function(model, id_states, delta_z= 0.5, dz_x=1,dz_y=1 )
{
  num_conformational_changes <- nrow(model$scheme$conformational_changes)
  num_domains <- nrow(model$scheme$domains)
  num_states <- length(model$states_and_index$conformational_states)


  stopifnot("more than 2 conformational changes, future releases will cover more changes" = num_conformational_changes <
              3)

  dcount <- model$scheme$domains %>% dplyr::count(id_conformational_change)

  domain_state <- purrr::reduce(id_states, function(x, i)
    c(
      x,
      model$states_and_index$conformational_states[[i]]$state_vector$changes_state
    ), .init = c())

  id_state <- rep(id_states, each = num_domains)

  dstates <- data.frame(id_state, model$scheme$domains, domain_state)

   n_states<-dstates %>% dplyr::group_by(id_state, id_conformational_change) %>%
    dplyr::summarise(n =sum(domain_state - 1)) %>%
    tidyr::pivot_wider(names_from =   id_conformational_change, values_from = n)%>%dplyr::arrange(`1`,`2`)%>%
     dplyr::ungroup()%>%dplyr::add_count(`1`,`2`)

   n_states$id_sub_state<-purrr::reduce(2:num_states,function(x,i) c(x, ifelse((n_states$`1`[i-1]==n_states$`1`[i])&&(n_states$`2`[i-1]==n_states$`2`[i]),x[length(x)]+1,1)),.init = c(1))->acc


   n_states$sub_pos= delta_z * (n_states$id_sub_state-(n_states$n + 1)/2) /n_states$n

   n_states$x=n_states$`1`+ dz_x* n_states$sub_pos
   n_states$y=n_states$`2`+ dz_y* n_states$sub_pos





  return (ggplot2::geom_point(data = n_states,ggplot2::aes(x=x, y=y, color=as.factor(id_state))))
}


