







#' Title
#'
#' @param conformational_changes
#' @param agonists
#'
#' @return
#' @export
#'
#' @examples
conformationalChanges <- function(conformational_changes,
                                  conformational_labels,
                                  agonists) {
  return (
    data.frame(
      id_conformational_change = 1:length(conformational_changes),
      conformational_change = conformational_changes,
      agonist = agonists
    )
  )
}

#' Title
#'
#' @param conformational_domain_changes
#'
#' @return
#' @export
#'
#' @examples
conformationaldomains <- function(conformational_domain_changes) {
  return (
    data.frame(
      id_domain = 1:length(conformational_domain_changes),
      id_conformational_change = conformational_domain_changes
    )
  )
}




#' Title
#'
#' @param labels
#' @param players_ids
#' @param domains_ids
#'
#' @return
#' @export
#'
#' @examples
conformationalInteractions <- function(labels, players_ids, domains_ids)
{
  interactions = data.frame(
    id_interaction = 1:length(labels),
    conformational_interaction = labels
  )

  interactions$players_ids = players_ids

  interactions$domains_ids = domains_ids


  return (interactions)
}

#' Title
#'
#' @param labels
#' @param players_ids
#' @param domains_ids
#'
#' @return
#' @export
#'
#' @examples
conductanceInteractions <- function(labels, players_ids, domains_ids)
{
  interactions = data.frame(
    id_interaction = 1:length(labels),
    conformational_interaction = labels
  )

  interactions$players_ids = players_ids

  interactions$domains_ids = domains_ids


  return (interactions)
}

#' Title
#'
#' @param id_conductance_interaction
#' @param kind
#' @param leakeage_label
#'
#' @return
#' @export
#'
#' @examples
conductanceInteractionInfo <- function(id_conductance_interaction,
                                       kind,
                                       leakeage_label)
{
  return (
    data.frame(
      id_conductance_interaction = id_conductance_interaction,
      kind = kind,
      leakeage_label = leakeage_label
    )
  )

}

#' Title
#'
#' @param id_conformational_change
#' @param id_interaction
#' @param pos_within_interaction
#' @param count
#'
#' @return
#' @export
#'
#' @examples
domainInteractionState <- function(id_conformational_change,
                                   id_interaction,
                                   pos_within_interaction,
                                   count)
{
  return (
    data.frame(
      id_conformational_change = id_conformational_change,
      id_interaction = id_interaction,
      pos_within_interaction = pos_within_interaction,
      count = count
    )
  )

}


#' Title
#'
#' @param id_conformational_change
#' @param id_interaction
#' @param pos_within_interaction
#' @param count
#'
#' @return
#' @export
#'
#' @examples
standardStates <- function(id_conformational_change,
                           id_interaction,
                           pos_within_interaction,
                           count)
{
  ss <- data.frame(
    id_conformational_change = id_conformational_change,
    id_interaction = id_interaction,
    pos_within_interaction = pos_within_interaction,
    count = count
  )
  return(ss)
}








#' Title
#'
#' @param domains
#' @param conformational_changes
#' @param conductaces
#' @param interactions
#'
#' @return
#' @export
#'
#' @examples
conformationalModelScheme <- function(conformational_changes,
                                      domains,
                                      interactions,
                                      conductaces)
{
  return (
    list(
      conformational_changes = conformational_changes,
      domains = domains,
      interactions = interactions,
      conductaces = conductaces
    )
  )
}




#' Title
#'
#' @param number_of_domains
#' @param expanded_state_number
#'
#' @return
#' @export
#'
#' @examples
conformationalChangeState <- function(scheme, expanded_state_number) {
  out = numeric(0)
  number_of_domains = nrow(scheme$domains)
  for (i in 0:number_of_domains - 1)
    out[i + 1] = (bitwAnd(expanded_state_number - 1, 2 ^ i) ==  (2 ^ i)) +
    1
  return (out)

}




#' Title
#'
#' @param Conformational_model_scheme
#' @param Conformational_change_state_vector
#'
#' @return
#' @export
#'
#' @examples
conformationalInteractionState <- function(scheme, state) {
  number_of_domains = nrow(scheme$domains)
  v_inter = scheme$interactions
  number_inter = nrow(v_inter)
  id_domain <- c()
  id_domain_2 <- c()
  id_conformational_change <- c()
  pos_within_interaction <- c()
  interaction_state <- c()

  id_conformational_change_2 <- c()
  id_interactions_2 <- c()
  pos_within_interaction_2 <- c()

  for (i in seq_len(number_of_domains)) {
    for (j in seq_len(number_inter)) {
      v_ipos = v_inter$domains_ids[[j]]

      for (k in seq_len(length(v_ipos))) {
        includes_i = FALSE
        includes_all_other = TRUE
        i_sub_position = 0
        for (kk in seq_len(length(v_ipos[[k]]))) {
          if (v_ipos[[k]][[kk]] == i) {
            includes_i = TRUE

            i_sub_position = kk

          } else {
            includes_all_other =
              includes_all_other && (state[v_ipos[[k]][[kk]]] == 2)
          }
        }
        if (includes_i)
        {
          id_domain <- append(id_domain, i)
          pos_within_interaction <-
            append(pos_within_interaction,
                   paste(j, i_sub_position, sep = "_"))
          interaction_state <- append(interaction_state, as.integer(includes_all_other))

          if (includes_all_other)
          {
            id_conformational_change_2 <- append(
              id_conformational_change_2,
              scheme$domains$id_conformational_change[i]
            )
            id_domain_2 <- append(id_domain_2, i)
            id_interactions_2 <- append(id_interactions_2, j)
            pos_within_interaction_2 <- append(pos_within_interaction_2, i_sub_position)

          }

        }
      }
    }
  }
  d <-
    data.frame(id_domain, pos_within_interaction, interaction_state)

  #dplyr::group_by(d, dplyr::across(dplyr::all_of(colnames(d)))) -> d
  #d <- dplyr::summarise(d, count = dplyr::n())
  d <-  tidyr::pivot_wider(
    d,
    names_from = c("pos_within_interaction"),
    values_from = c("interaction_state"),
    values_fill = 0
  )
  d2 <- data.frame()

  if (length(id_interactions_2) > 0)
  {
    d2 <-
      data.frame(
        id_domain = id_domain_2,
        #  id_conformational_change,
        id_interaction = id_interactions_2,
        pos_within_interaction = pos_within_interaction_2
      )
    dplyr::group_by(d2, dplyr::across(dplyr::all_of(colnames(d2)))) -> d2
    d2 <- dplyr::summarise(d2, count = dplyr::n())
  }
  else
  {
    d2 <- data.frame(matrix(nrow = 0, ncol = 4))
    colnames(d2) <- c("id_domain",
                      "id_interaction",
                      "pos_within_interaction",
                      "count")
  }
  return (list(matrix_form = d, list_form = d2))
}


conformationalInteractionState_old <- function(scheme, state) {
  number_of_domains = nrow(scheme$domains)
  v_inter = scheme$interactions
  number_inter = nrow(v_inter)
  id_domain <- c()
  id_conformational_change <- c()
  id_interactions <- c()
  pos_within_interaction <- list()
  for (i in seq_len(number_of_domains)) {
    for (j in seq_len(number_inter)) {
      v_ipos = v_inter$domains_ids[[j]]

      for (k in seq_len(length(v_ipos))) {
        includes_i = FALSE
        includes_all_other = TRUE
        i_sub_position = 0
        for (kk in seq_len(length(v_ipos[[k]]))) {
          if (v_ipos[[k]][[kk]] == i) {
            includes_i = TRUE

            i_sub_position = kk

          } else {
            includes_all_other =
              includes_all_other && (state[v_ipos[[k]][[kk]]] == 2)
          }
        }
        if (includes_all_other && includes_i)
        {
          id_conformational_change <- append(id_conformational_change,
                                             scheme$domains$id_conformational_change[i])
          id_domain <- append(id_domain, i)
          id_interactions <- append(id_interactions, j)
          pos_within_interaction <- append(pos_within_interaction, i_sub_position)

        }
      }
    }
  }
  d <-
    data.frame(id_domain, #  id_conformational_change,
               id_interaction = id_interactions, pos_within_interaction)
  dplyr::group_by(d, dplyr::across(dplyr::all_of(colnames(d)))) -> d
  d <- dplyr::summarise(d, count = dplyr::n())
  # d <-  tidyr::pivot_wider(
  #   d,
  #   names_from = c("id_interactions", "pos_within_interaction"),
  #   values_from = c("count"),
  #   values_fill = 0
  # )
  return (d)
}


#' Title
#'
#' @param scheme
#' @param expanded_state_number
#' @param change
#'
#' @return
#' @export
#'
#' @examples
conformationalState <- function(scheme,
                                expanded_state_number = NULL,
                                change = NULL) {
  if (is.null(change))
    change = conformationalChangeState(scheme = scheme, expanded_state_number = expanded_state_number)

  inter = conformationalInteractionState(scheme, change)

  return (list(changes_state = change, interactions_state = inter))
}







#' Title
#'
#' @param const
#' @param const
#'
#' @return
#' @export
#'
#' @examples
to_state_count <- function(scheme, state) {
  stopifnot (length(state$changes_state)
             == nrow(scheme$domains))
  v_state = state$changes_state
  d = state$interactions_state$matrix_form
  d$conformational_state = v_state
  d <- dplyr::group_by(d, dplyr::across(-id_domain))
  d <- dplyr::summarise(d, count = dplyr::n())
  return(d)
}





#' Title
#'
#' @param scheme
#' @param state
#'
#' @return
#' @export
#'
#' @examples
to_state_conductance_count <- function(scheme, state) {
  number_of_domains = nrow(scheme$domains)
  v_inter = scheme$conductaces
  number_inter = nrow(v_inter)
  id_conductance_interaction <- c()
  conductance_interaction_state <- c()
  for (j in seq_len(number_inter)) {
    v_ipos = v_inter$domains_ids[[j]]
    for (k in seq_len(length(v_ipos))) {
      includes_all_other = TRUE
      for (kk in seq_len(length(v_ipos[[k]]))) {
        includes_all_other =
          includes_all_other && (state[v_ipos[[k]][[kk]]] == 2)
      }

      id_conductance_interaction <- append(id_conductance_interaction, j)
      conductance_interaction_state <- append(conductance_interaction_state,
                                              as.integer(includes_all_other))

    }
  }
  d <-
    data.frame(id_conductance_interaction, conductance_interaction_state)

  dplyr::group_by(d, id_conductance_interaction) -> d
  d <- dplyr::summarise(d, count = sum(conductance_interaction_state))
  # d <-  tidyr::pivot_wider(
  #   d,
  #   names_from = c("pos_within_interaction"),
  #   values_from = c("interaction_state"),
  #   values_fill = 0
  # )
  return (d)
}









to_state_conductance_count_old <- function(scheme, state) {
  v_state = state$changes_state

  v_inter = scheme$conductaces
  number_inter = nrow(v_inter)

  id_conductance_interaction = c()
  active = c()
  for (i in seq_len(number_inter)) {
    pos_domains = v_inter$domains_ids[[i]]
    players_id = v_inter$players_ids[[i]]
    for (k  in seq_len(length(pos_domains)) ){
      includes_all = TRUE
      for (kk in seq_len(length(pos_domains[[k]]))) {
        includes_all = includes_all && (v_state[players_id[kk]] == 2)
      }
      if (includes_all)

        id_conductance_interaction <- append(id_conductance_interaction, i)

    }
  }
  d <- data.frame(id_conductance_interaction)
  d <- dplyr::group_by(d, id_conductance_interaction)
  d <- dplyr::summarise(d, count = dplyr::n())



  return (d)

}





#' Title
#'
#' @param state_count
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#' @examples
to_state_barcode <-  function(state_count)
{
  r_state_barcode <- state_count %>% tidyr::pivot_wider(names_from = tidyr::everything(),
                                                        values_from = count) %>%
    colnames() %>% paste(collapse = "__")

}



#' Title
#'
#' @param scheme
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_wider
#' @return
#' @export
#'
#' @examples
conformational_states_and_index <- function (scheme) {
  number_of_domains = nrow(scheme$domains)


  id_state = c()
  id_state_expanded = c()
  conformational_states = list()
  state_barcode = c()

  for (i in seq_len(2 ^ number_of_domains)) {
    state = conformationalState(scheme = scheme, expanded_state_number = i)

    state_count = to_state_count(scheme, state)

    state_conductance =
      to_state_conductance_count(scheme, state$changes_state)



    r_state_barcode <- to_state_barcode(state_count)

    if (!r_state_barcode %in% state_barcode)
    {
      state_barcode <- append(state_barcode, r_state_barcode)
      i_state = length(state_barcode)
      id_state <- append(id_state, i_state)
      id_state_expanded <- append(id_state_expanded, i)
      conformational_states[[i_state]] <-
        list(
          state_count = state_count,
          stabarcode = r_state_barcode,
          state_vector = state,
          state_conductance = state_conductance
        )


    }


  }
  barcode_map = id_state
  names(barcode_map) <- state_barcode

  return (
    list(
      id_state = id_state,
      id_state_expanded = id_state_expanded,
      barcode_map = barcode_map,
      conformational_states = conformational_states
    )
  )
}




#' Title
#'
#' @param state
#' @param ith_domain
#'
#' @return
#' @export
#'
#' @examples
change_state_conformation <- function(state, ith_domain) {
  out <- state

  if (out[ith_domain] == 1)
    out[ith_domain] = 2
  else
    out[ith_domain] = 1

  return (out)

}


#' Title
#'
#' @param scheme
#' @param state
#' @param ith_domain
#'
#' @return
#' @export
#'
#' @examples
change_conformation <- function(scheme, state, ith_domain) {
  v_change = change_state_conformation(state = state, ith_domain = ith_domain)
  return (conformationalState(scheme = scheme, change = v_change))
}




#' Title
#'
#' @param states_and_index
#' @param j_count
#'
#' @return
#' @export
#'
#' @examples
to_state_index <- function(states_and_index, j_count)
{
  j_barcode <- to_state_barcode(j_count)
  return (states_and_index$barcode_map[j_barcode][[1]])
}


#' Title
#'
#' @param scheme
#' @param states_and_index
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
conformational_transition_list <- function(scheme, states_and_index) {
  #std::vector<std::vector<Conformational_transition>> out(states().size());

  number_of_states = length(states_and_index$id_state)
  #auto scheme = get<Conformational_change_scheme>(model());

  states = states_and_index$conformational_states

  out <- list()
  for (i in seq_len(number_of_states)) {
    #//  auto v_count=get<Conformational_state_count>(states()[i]);
    #v_state = get<Conformational_state_vector>(states()[i]);

    v_state = states[[i]]$state_vector

    # auto v_change = get<Conformational_change_state_vector>(v_state());

    v_change_state = v_state$changes_state

    #auto v_interactions =
    #  get<Conformational_interactions_state_vector>(v_state());

    v_interactions = v_state$interactions_state$list_form

    i_start = i


    tran <- list()

    #std::vector<Conformational_transition> tran;
    for (j in  1:length(v_change_state)) {
      j_change = scheme$domains$id_conformational_change[j]

      #    j_inter = v_interactions %>% filter(id_domain == j)

      j_state = change_conformation(scheme, v_change_state, j)

      #       auto Maybe_j_count = to_state_count(model, j_state);
      #       if (!Maybe_j_count)
      #         return Maybe_j_count.error();
      #       else {
      #         auto &j_count = Maybe_j_count.value();
      j_count <- to_state_count(scheme, j_state)


      # auto &j_count = Maybe_j_count.value();
      # auto Maybe_j_index = to_state_index(map, j_count);
      # if (!Maybe_j_index)
      #   return Maybe_j_index.error();
      # else {
      #   auto i_end = Conformational_transition_landing_state_index(
      #     Maybe_j_index.value());

      j_index = to_state_index(states_and_index = states_and_index, j_count)


      if (!j_index %in% names(tran)) {
        #     auto a = get<Agonist_dependency>(j_change());
        #     auto d = Conformational_transition_direction(!v_change()[j]());
        #     auto l = get<Conformational_change_label>(j_change());
        #     auto m = Conformational_transition_mulitplicity(1);
        #     tran.push_back(Conformational_transition(
        #       Vector_Space(i_start, i_end, a, d, l, j_inter, m)));
        #   }
        #   #           }
        tran[[paste0(j_index)]] <- list(
          i_start = i,
          i_end = j_index,
          agonist_dependency = scheme$conformational_changes$agonist[j_change],
          conformational_transition_direction = j_state$changes_state[j] -
            v_change_state[j],
          id_conformational_change = scheme$conformational_changes$id_conformational_change[j_change],
          conformational_change = scheme$conformational_changes$conformational_change[j_change],
          conformational_interactions = dplyr::filter(j_state$interactions_state$list_form, id_domain ==
                                                        j),
          conformational_transition_mulitplicity = 1
        )
      }
      else {
        tran[[paste0(j_index)]]$conformational_transition_mulitplicity =
          tran[[paste0(j_index)]]$conformational_transition_mulitplicity +
          1
      }


    }
    out[[i]] <- tran
  }

  return(out)


}








#' Title
#'
#' @param conformational_changes
#' @param domains
#' @param interactions
#' @param conductances
#' @param conductance_info
#' @param standard_states
#'
#' @return
#' @export
#'
#' @examples
ConformationalModel <- function(conformational_changes,
                                domains,
                                interactions,
                                conductances,
                                standard_states,
                                conductance_info)
{
  scheme = conformationalModelScheme(conformational_changes, domains, interactions, conductances)
  states_and_index = conformational_states_and_index(scheme)

  transitions = conformational_transition_list(scheme, states_and_index)

  return (
    list(
      scheme = scheme,
      states_and_index = states_and_index,
      transitions = transitions,
      standard_states = standard_states,
      conductance_info = conductance_info
    )
  )
}




#' Title
#'
#' @param inter
#' @param par
#' @param trr
#' @param v_int
#'
#' @importFrom magrittr %>%
#' @return
#' @export
#'
#' @examples
calc_Qij <- function (interactions,
                      parameters,
                      trr,
                      standard_states)
{
  chla <- trr$conformational_change
  id_chla <- trr$id_conformational_change
  conformational_interactions <- trr$conformational_interactions %>% dplyr::ungroup() %>%
    dplyr::select (!id_domain)
  st <- dplyr::select(standard_states, !id_conformational_change)
  st$count = -st$count
  conformational_interactions <- dplyr::bind_rows(conformational_interactions, st) %>%
    dplyr::group_by(id_interaction, pos_within_interaction) %>% dplyr::summarise(count =
                                                                                   sum(count)) %>% dplyr::filter(count != 0)


  #   get<Conformational_interactions_domain_state>(tr()) - st[chla];


  ag = trr$agonist_dependency

  #auto d = get<Conformational_transition_direction>(tr());
  #auto chla = get<Conformational_change_label>(tr());
  #auto n = get<Conformational_transition_mulitplicity>(tr());


  d <- trr$conformational_transition_direction




  n <- trr$conformational_transition_mulitplicity



  # auto Maybe_i_base = d() ? names[chla() + "_on"] : names[chla() + "_off"];
  #
  # if (!Maybe_i_base)
  #   return Maybe_i_base.error();
  # auto i_base = Maybe_i_base.value();
  base <- ifelse(d > 0, parameters[paste0(chla, "_on")], parameters[paste0(chla, "_off")])

  # auto out = n() * par()[i_base];
  out = n * base

  # for (auto ii = v_int().begin(); ii != v_int().end(); ++ii) {
  for (ii in seq_len(nrow(conformational_interactions))) {
    #   auto factor_la = get<Conformational_interaction_label>(
    #     inter()[get<Conformational_interaction_index>(ii->first)()]())();
    factor_la = interactions$conformational_interaction[conformational_interactions$id_interaction[ii]]
    factor_ipos =  conformational_interactions$pos_within_interaction[ii]
    factor_power = conformational_interactions$count[ii]


    #   auto Maybe_i_Factor = names[factor_la];

    Factor <- parameters[factor_la]


    #   auto Maybe_i_Factor_pos =
    #     names[factor_la + "_" + std::to_string(factor_ipos)];
    Factor_pos <- parameters[paste0(factor_la, "_", factor_ipos-1)]


    #   if (!Maybe_i_Factor || !Maybe_i_Factor_pos)
    #     return error_message(Maybe_i_Factor.error()() +
    #                            Maybe_i_Factor_pos.error()());
    #   else {
    #     auto i_Factor = Maybe_i_Factor.value();
    #     auto i_Factor_pos = Maybe_i_Factor_pos.value();
    #
    #     using std::pow;
    #     if (d())
    #       out = out * pow(par()[i_Factor_pos], factor_power);
    #     else
    #       out = out * pow(par()[i_Factor_pos] / par()[i_Factor], factor_power);

    out = ifelse(d > 0,
                 out * Factor_pos ^ factor_power,
                 out * Factor_pos / (Factor ^ factor_power))

  }
  return (out)
}














#' Title
#'
#' @param model
#' @param names
#' @param par
#'
#' @return
#' @export
#'
#' @examples
make_Q0_Qa <- function(model, parameters) {
  N <- length(model$states_and_index$id_state)
  inter <- model$scheme$interactions

  tr <- model$transitions
  stopifnot(length(tr) == N)

  v_Q0 <- matrix(data = 0,
                 nrow = N,
                 ncol = N)
  v_Qa <- matrix(data = 0,
                 nrow = N,
                 ncol = N)

  for (i  in seq_len(N)) {
    for (j in seq_len(length(tr[[i]]))) {
      #  auto trr = tr()[i][j];

      trr = tr[[i]][[j]]

      #auto v_i = get<Conformational_transition_initiating_state_index>(trr());
      v_i = trr$i_start
      stopifnot(v_i == i)
      v_j = trr$i_end
      ag = trr$agonist_dependency
      d = trr$conformational_transition_direction

      # Maybe_error<var::Op_t<transformation_type_t<P>, double>> Maybe_qij;
      # if constexpr (std::is_same_v<Conformational_model_,
      #               Conformational_model_standarized>)
      # Maybe_qij =
      #   calc_Qij<Id>(inter, names, par, trr,
      #                get<Conformation_change_standard_map>(model()));
      # else
      #   Maybe_qij = calc_Qij<Id>(inter, names, par, trr);

      v_qij = calc_Qij(
        scheme$interactions,
        parameters = parameters,
        trr = trr,
        standard_states = standard_states
      )

      if (trr$agonist_dependency == "" || d < 0) {
        #   // set(v_Q0(), i, i, v_Q0()(i, i) - Maybe_qij.value());  later change it
        v_Q0[i, v_j] <- v_qij
        #   // back
        #   set(v_Q0(), i, v_j()(), std::move(Maybe_qij.value()));
      } else {
        #   //  set(v_Qa(), i, i, v_Qa()(i, i) - Maybe_qij.value());  same
        #   set(v_Qa(), i, v_j()(), std::move(Maybe_qij.value()));
        v_Qa[i, v_j] <- v_qij

      }
    }
  }
  return (list(Q0 = v_Q0, Qa = v_Qa))
}









make_model <- function (model, parameters)
{
  Q0Qa = make_Q0_Qa(model, parameters)

  # g = impl::make_g<Id>(model, names, p);
  # if (!Maybe_Q0Qa || !Maybe_g)
  #   return error_message(Maybe_Q0Qa.error()() + Maybe_g.error()());
  # else {
  #   auto [v_Q0, v_Qa] = std::move(Maybe_Q0Qa.value());
  #   auto v_g = std::move(Maybe_g.value());
  #
  #   return std::tuple(std::move(v_Q0), std::move(v_Qa), std::move(v_g));
  # }
}
