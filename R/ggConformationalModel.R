










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
  angles = (seq_len(n_sides) - 1) * 2 * pi / n_sides
  xi = radius * sin(angles) * h_x
  yi =  radius * cos(angles) * h_y

  xs = x + xi * cos(rotation) - yi * sin(rotation)
  ys =  y + xi * sin(rotation) + yi * cos(rotation)
  return (data.frame(x = xs, y = ys))
}


half_arrow <- function(x0,
                       y0,
                       xend,
                       yend,
                       head_size,
                       head_ange,
                       gap,
                       shift)
{
  dx = x0 - xend
  dy = y0 - yend

  total_length = sqrt(dx ^ 2 + dy ^ 2)

  rx = dx / total_length
  ry = dy / total_length


  stopifnot(total_length > 2 * gap + head_size)

  x0c = x0 - gap * rx + shift * ry
  y0c = y0 - gap * ry - shift * rx
  xendc = xend + gap * rx + shift * ry
  yendc = yend + gap * ry - shift * rx


  head_edge_length = head_size / cos(head_ange)

  xa = (rx * cos(head_ange) + ry * sin(head_ange))  * head_edge_length +
    xendc
  ya = (-rx * sin(head_ange) + ry * cos(head_ange)) * head_edge_length +
    yendc

  return (data.frame(x = c(x0c, xendc, xa), y = c(y0c, yendc, ya)))
}






#' Title
#'
#' @param Q
#'
#' @return
#' @export
#'
#' @examples
Q_to_dataframe <- function(Q)
{
  n = nrow(Q)
  i_m <- matrix(rep(1:n, n), nrow = n)
  j_m <- matrix(rep(1:n, each = n), nrow = n)

  return(data.frame(i = i_m[Q > 0], j = j_m[Q > 0], qij = Q[Q > 0]))

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
    fill_i [ind] <- fill [min(length(fill), i)]
  }


  d <- data.frame(
    x = xi,
    y = yi,
    gi = gi,
    color = color_i,
    fill = fill_i
  )

  return (ggplot2::geom_polygon(
    data = d,
    ggplot2::aes(
      x = x,
      y = y,
      group = gi,
      fill = fill,
      color = color
    )
  ))
}





#' Title
#'
#' @param x
#' @param y
#' @param xend
#' @param yend
#' @param head_size
#' @param head_angle
#' @param alpha
#' @param linetype
#' @param linewidth
#' @param color
#' @param gap
#' @param shift
#'
#' @return
#' @export
#'
#' @examples
geom_half_arrow <- function(x,
                            y,
                            xend,
                            yend,
                            head_size,
                            head_angle,
                            gap,
                            shift,
                            color,
                            alpha,
                            linetype,
                            linewidth) {
  stopifnot(length(x) == length(y))
  stopifnot(length(xend) == length(y))
  stopifnot(length(yend) == length(y))
  stopifnot(length(color) == 1 | length(color) == length(x))
  stopifnot(length(alpha) == 1 | length(alpha) == length(x))
  stopifnot(length(linetype) == 1 | length(linetype) == length(x))
  stopifnot(length(linewidth) == 1 | length(linewidth) == length(x))

  n = 3 * length(x)

  xi = vector("numeric", n)
  yi = vector("numeric", n)
  gi = vector("numeric", n)
  color_i = ifelse(length(color)>1,rep(color[1], n), color)
  alpha_i = ifelse(length(alpha)>1,rep(alpha[1], n), alpha)
  linetype_i = ifelse(length(linetype)>1,rep(linetype[1], n), linetype)
  linewidth_i = ifelse(length(linewidth)>1,rep(linewidth[1], n), linewidth)



  r_i = 0
  for (i in seq_len(length(x))) {
    ind <- 1:3 + (i - 1) * 3
    v <- half_arrow(
      x0 = x[i],
      y0 = y[i],
      xend = xend[i],
      yend = yend[i],
      head_size = head_size,
      head_ange = head_angle,
      gap = gap,
      shift = shift
    )
    xi[ind] <- v$x
    yi[ind] <- v$y
    gi[ind] <-  i
    if (length(color)>1) color_i[ind] <- color [min(length(color), i)]
    if (length(alpha)>1) alpha_i[ind] <- alpha [min(length(alpha), i)]
    if (length(linetype)>1) linetype_i[ind] <- linetype [min(length(linetype), i)]
    if (length(linewidth)>1) linewidth_i[ind] <- linewidth [min(length(linewidth), i)]

  }


  d <- data.frame(
    x = xi,
    y = yi,
    group = gi,
    color = color_i,
    alpha = alpha_i,
    linetype = linetype_i,
    linewidth = linewidth_i
  )



  e<-expression(ggplot2::geom_path(
    data = d,
    ggplot2::aes(
      x = x,
      y = y,
      group = group,
    )
  ))

  if (length(color)>1)
    e[[1]][[3]]$color=color_i
  else
    e[[1]]$color=color

  if (length(alpha)>1)
    e[[1]][[3]]$alpha=alpha_i
  else
    e[[1]]$alpha=alpha

  if (length(linetype)>1)
    e[[1]][[3]]$linetype=linetype_i
  else
    e[[1]]$linetype=linetype

  if (length(linewidth)>1)
    e[[1]][[3]]$linewidth=linewidth_i
  else
    e[[1]]$linewidth=linewidth




  return (eval(e))
}












#' Title
#'
#' @param model
#' @param id_states
#' @param delta_z
#' @param dz_x
#' @param dz_y
#' @param state_radius
#'
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
g_ConformationalStates_data <- function(model,
                                        id_states,
                                        state_radius = 0.05,
                                        delta_z = 0.5,
                                        dz_x = 1,
                                        dz_y = 1)
{
  n_states <- model$domain_states %>% dplyr::filter(id_state %in% id_states)
  num_states <- length(id_states)

  num_domains <- nrow(model$scheme$domains)

  domain_angle <- 0:(num_domains - 1) / num_domains * 2 * pi

  domain_x <- cos(domain_angle) * state_radius
  domain_y <- sin(domain_angle) * state_radius




  n_states$sub_pos = delta_z * (n_states$id_sub_state - (n_states$n + 1) /
                                  2) / n_states$n

  conf_1 = model$scheme$conformational_changes$conformational_change[1]
  conf_2 = model$scheme$conformational_changes$conformational_change[2]


  n_states$x = n_states[[conf_1]] + dz_x * n_states$sub_pos + domain_x[n_states$id_domain]
  n_states$y = n_states[[conf_2]] + dz_y * n_states$sub_pos + domain_y[n_states$id_domain]

  n_states$n_sides = ifelse(n_states$agonist == "",
                            64,
                            ifelse(n_states$domain_state == 1, 1, 3))
  n_states$radius = ifelse(n_states$agonist == "", state_radius * 0.8, state_radius *
                             0.8)

  n_states$rotation = rep(domain_angle - pi / 6, num_states)

  n_states$fill = ifelse(n_states$domain_state == 1, "off", "on")

  return (n_states)

}





#' Title
#'
#' @param delta_z
#' @param dz_x
#' @param dz_y
#' @param n_states
#' @param state_radius
#'
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
g_ConformationalStates <- function(n_states,
                                   state_radius = 0.05,
                                   delta_z = 0.5,
                                   dz_x = 1,
                                   dz_y = 1)
{
  return (
    ggplot2::ggplot() + geom_regular_polygon(
      x = n_states$x,
      y = n_states$y,
      n_sides = n_states$n_sides,
      rotation = n_states$rotation,
      radius = n_states$radius,
      color = 1,
      fill = n_states$fill,
      h_x = 1,
      h_y = 1
    ) +
      ggplot2::scale_color_manual(
        aesthetics = c("fill"),
        values = c("white", "black")
      ) + ggplot2::coord_fixed()
  )
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
g_Rate_arrows <- function(model,
                          n_domains,
                          Q0,
                          Qa,
                          id_states,
                          state_radius = 0.05,
                          delta_z = 0.5,
                          dz_x = 1,
                          dz_y = 1)
{
  conformational_changes_list <- model$scheme$conformational_changes$conformational_change

  dQ0 <- Q_to_dataframe(Q0[id_states,id_states])

  dQa <- Q_to_dataframe(Qa[id_states,id_states])

  dQ <- dplyr::bind_rows(dQ0, dQa)

  n_states <- n_domains %>% dplyr::filter(id_domain == 1) %>%
    dplyr::select(id_state,
                  dplyr::all_of(conformational_changes_list),
                  sub_pos)

  conf_1 = conformational_changes_list[1]
  conf_2 = conformational_changes_list[2]


  n_states$x = n_states[[conf_1]] + dz_x * n_states$sub_pos
  n_states$y = n_states[[conf_2]] + dz_y * n_states$sub_pos

  n_states <- n_states %>% dplyr::select(id_state, x, y)


  dQ <- dplyr::left_join(
    dplyr::left_join(dQ, n_states, by = dplyr::join_by(i == id_state)),
    n_states,
    by = dplyr::join_by(j == id_state),
    suffix = c(".i", ".j")
  )




  return (
    geom_half_arrow(
      x = dQ$x.i,
      y = dQ$y.i,
      xend = dQ$x.j,
      yend = dQ$y.j,
      head_size = state_radius * 1,
      head_angle = pi / 2 / 3,
      gap = state_radius * 2,
      shift = state_radius * 0.2,
      color = dQ$qij,
      alpha = 1,
      linetype = 1,
      linewidth = 1
    )
  )
}


#' Title
#'
#' @param model
#' @param Q0
#' @param Qa
#' @param id_states
#' @param state_radius
#' @param delta_z
#' @param dz_x
#' @param dz_y
#'
#' @return
#' @export
#'
#' @examples
g_States_and_Rates <- function(model,
                               Q0,
                               Qa,
                               id_states,
                               state_radius = 0.05,
                               delta_z = 0.5,
                               dz_x = 1,
                               dz_y = 1)
{
  n_domains <- g_ConformationalStates_data(
    model,
    id_states = id_states,
    state_radius = state_radius,
    delta_z = delta_z,
    dz_x = 1,
    dz_y = 1
  )

  return(
    g_ConformationalStates(
      n_states = n_domains,
      state_radius = state_radius,
      delta_z = delta_z,
      dz_x = dz_x,
      dz_y = dz_y
    ) +
      g_Rate_arrows(model = model,
                    n_domains = n_domains,
                    Q0 = Q0,Qa = Qa,
                    id_states = id_states,
                    state_radius = state_radius,delta_z = delta_z,dz_x = dz_x,dz_y = dz_y))
    }
