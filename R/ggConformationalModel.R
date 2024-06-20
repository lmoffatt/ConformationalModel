






unit_circle <- function(n_steps)
{
  return (list(x, y))
}



regular_polygon <- function(x,
                            y,
                            n_sides,
                            radius,
                            rotation,
                            h_x = 1,
                            h_y = 1)
{
  angles = (seq_len(n_sides)-1) * 2 * pi / n_sides
  xi = radius * sin(angles) * h_x?
  yi =  radius * cos(angles) * h_y

  xs = x + xi * cos(rotation) - yi * sin(rotation)
  ys =  y + xi * sin(rotation) + yi * cos(rotation)
  return (data.frame(x = xs, y = ys))
}






#' Title
#'
#' @param x
#' @param y
#' @param radius
#' @param fill
#' @param color
#' @param n_steps
#'
#' @return
#' @export
#'
#' @examples
geom_regular_polygon <- function(x,
                                 y,
                                 n_sides,
                                 radius,
                                 rotation,
                                 h_x,
                                 h_y,
                                 color,
                                 fill) {
  stopifnot(length(x) == length(y))
  stopifnot(length(radius) == 1 | length(radius) == length(x))
  stopifnot(length(n_sides) == 1 | length(n_sides) == length(x))
  stopifnot(length(rotation) == 1 | length(rotation) == length(x))
  stopifnot(length(h_x) == 1 | length(h_x) == length(x))
  stopifnot(length(h_y) == 1 | length(h_y) == length(x))

  stopifnot(length(color) == 1 | length(color) == length(x))
  stopifnot(length(fill) == 1 | length(fill) == length(x))


  n = ifelse(length(n_sides) == 1, n_sides * length(x), sum(n_sides))

  xi = vector("numeric", n)
  yi = vector("numeric", n)
  gi = vector("numeric", n)
  color_i = rep(color[1], n)
  fill_i = rep(fill[1], n)



  r_i = 0
  for (i in seq_len(length(x))) {
    nsides = n_sides [min(length(n_sides), i)]
    ind <- 1:nsides + r_i
    r_i <- r_i + nsides
    v <- regular_polygon(
      x[i],
      y[i],
      n_sides = nsides,
      radius = radius [min(length(radius), i)],
      rotation = rotation [min(length(rotation), i)],
      h_x = h_x [min(length(h_x), i)],
      h_y = h_y [min(length(h_y), i)]
    )
    xi[ind] <- v$x
    yi[ind] <- v$y
    gi[ind] <-  i
    color_i[ind] <- color [min(length(color), i)]
    fill_i [ind]<- fill [min(length(fill), i)]
  }


  d <- data.frame(
    x = xi,
    y = yi,
    gi = gi,
    color = color_i,
    fill = fill_i
  )

  return (ggplot2::geom_polygon(data = d, ggplot2::aes(
    x = x,
    y = y,
    group = gi,
    fill = fill,
    color = color
  )))
}














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
g_ConformationalStates <- function(model, id_states,state_radius=0.05, delta_z= 0.5, dz_x=1, dz_y=1 )
{
  n_states <- model$domain_states%>%dplyr::filter(id_state %in% id_states)
  num_states<-length(id_states)

  num_domains<-nrow(model$scheme$domains)

  domain_angle<- 0:(num_domains-1) /num_domains *2 * pi

  domain_x<- cos(domain_angle)* state_radius
  domain_y<- sin(domain_angle)* state_radius




  n_states$sub_pos= delta_z * (n_states$id_sub_state-(n_states$n + 1)/2) /n_states$n

  conf_1=model$scheme$conformational_changes$conformational_change[1]
  conf_2=model$scheme$conformational_changes$conformational_change[2]


  n_states$x=n_states[[conf_1]]+ dz_x* n_states$sub_pos + domain_x[n_states$id_domain]
  n_states$y=n_states[[conf_2]]+ dz_y* n_states$sub_pos+ domain_y[n_states$id_domain]

  n_states$n_sides=ifelse(n_states$agonist=="",64,ifelse(n_states$domain_state==1,1,3))
  n_states$radius=ifelse(n_states$agonist=="",state_radius*0.8,state_radius*0.8)

  n_states$rotation=rep(domain_angle - pi/6, num_states)

  n_states$fill=ifelse(n_states$domain_state==1,"off","on")


  return (ggplot2::ggplot()+geom_regular_polygon(x = n_states$x, y=n_states$y, n_sides =n_states$n_sides,rotation = n_states$rotation,
                               radius= n_states$radius, color=1,fill=n_states$fill,h_x = 1,h_y = 1)+
            ggplot2::scale_color_manual(aesthetics = c( "fill"), values = c("white", "black"))+ggplot2::coord_fixed())
}


