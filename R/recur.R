#' Recurrent Survival Data Format
#'
#' @description Reformats recurrent event data (wide) into different models for
#'   survival analysis, but can also be used for simple survival analysis tables
#'   as well. The general format is how data tends to be collected. There is
#'   left and right censoring date, a labeled event column that contains the
#'   date of the event, and a censoring column for a final censoring event. The
#'   accepted parameter options are listed, with the type of table that will be
#'   generated:
#'
#'   * __traditional__: Traditional survival table that has single censoring
#'   event (`trad`)
#'
#'   * __counting__ : Formally called the Andersen and Gill model (`ag`).
#'   Counting process model assumes each event is independent and that a subject
#'   contributes to the risk set during the time under observation. Multiple
#'   events are treated as a new (but delayed) entry that is followed until the
#'   next event. This means subjects under observation are at risk for a second
#'   event, even without having had a prior event. There are thus no _strata_ in
#'   the model.
#'
#'   * __marginal__: Formally called the Wei-Lin-Weissfield model, but more
#'   commonly known as a marginal model (`marginal`). Marginal models assumes
#'   each event is a separate process. Each subject is at risk for all events.
#'   The time for an event starts at the beginning of follow-up for each
#'   subject. Thus, each risk period is considered a different _strata_
#'   (regardless of if subject had an event or not).
#'
#'   * __conditional A__: Formally called the Prentice, Williams, and Peterson
#'   total time model (`pwptt`). Conditional A models order events by
#'   stratification, based on the number of events prior. All subjects are at
#'   risk for the left _strata_, but only those with a previous event are at
#'   risk for a successive event. The total time to event is used.
#'
#'   * __conditional B__: Formally called the Prentice, Williams, and Peterson
#'   gap time model (`pwpgt`). Conditional B models also order events by strata
#'   (like conditional A), however the time to outcome is defined as the gap
#'   between the time of previous event.
#'
#' @details This function takes every event date, and creates several types of
#'   recurrent event tables. It orders the data chronologically for repeat
#'   events. Currently does normal (left event) and recurrent models (counting,
#'   marginal, and conditional A and B models). Further details can be found at
#'   [IDRE](https://stats.idre.ucla.edu/sas/faq/how-can-i-model-repeated-events-survival-analysis-in-proc-phreg/).
#'
#'   * For recurrent events, the final censoring event can include death, or can
#'   be ignored if its not considered a failure event.
#'
#'   * For traditional survival analysis, `censor` is required and `event_dates`
#'   should be left as NULL. The function will do the rest.
#'
#'   __Performance__: Importantly, for large datasets of recurrent data (>500
#'   rows), this function will show significant slow-down since it uses an
#'   intuitive approach on defining the datasets. Future iterations will create
#'   a vectorized approach that should provide performance speed-ups.
#'
#' @return A data frame organized into a survival table format. Output options
#'   are in __Details__. Generally, the following columns are generated:
#'
#'   - __id__: An ID column is created
#'
#'   - __start__: A formatted start time, usually 0
#'
#'   - __stop__: A formatted stop time, in days, from prior event
#'
#'   - __status__: If event occurred or not
#'
#'   - __strata__: Event strata that is being applied
#'
#' @param data A dataframe containing the subsequent parameters
#'
#' @param model_type Model type that is indicated:
#'
#'   * `trad` makes traditional survival table
#'
#'   * `ag` makes table with risk periods starting at time of prior event
#'   without conditional strata
#'
#'   * `marginal` makes table with risk periods from entry to censorship with
#'   strata per each event
#'
#'   * `pwptt` makes table with risk periods starting at time of prior event
#'   with conditional strata
#'
#'   * `pwpgt` makes table with risk periods of each time interval between
#'   events,  with conditional strata
#'
#' @param id Column in dataframe that contains unique IDs for each row
#'
#' @param left Column with left/enrollment dates
#'
#' @param right Column with right/censoring time point, or right contact
#'
#' @param censor Column that names if death/final censorship is known (0 or 1).
#'   The default is that, if no censorship information is given, that are no
#'   failure events at time of right contact. `censor` is not required for
#'   recurrent event analysis, but is required for traditional survival tables.
#'
#' @param event_dates Vector of columns that contain event dates
#'
#' @examples
#' \donttest{
#' # Data
#' data("mims")
#'
#' # Parameters
#' id <- "patid"
#' left <- "left_visit_date_bl"
#' right <- "ldka"
#' event_dates <- c("mi_date_1", "mi_date_2", "mi_date_3")
#' model_type <- "marginal"
#' censor <- "DEATH_CV_YN"
#'
#' # Run analysis
#' out <- recur(
#'   mims, model_type, id, left, right, censor, event_dates
#' )
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select all_of mutate arrange group_by ungroup
#' @importFrom purrr map
#' @importFrom tidyr nest unnest pivot_longer pivot_wider
#' @export
#' @rdname recur
recur <- function(data,
									model_type,
									id,
									left,
									right,
									censor = NULL,
									event_dates = NULL) {

	# Check to see if censoring column is available
	if (is.null(censor)) {
		data$censor <- 0
		warning("Censorship data was not provided.")
	}

	# Check on event dates
	if (model_type == "trad") {
		if (!is.null(event_dates)) {
			warning("Event dates shouldn't be specified for traditional survival table.")
		}
	} else if (model_type != "trad") {
		if (is.null(event_dates)) {
			stop("Event dates should be specified for recurrent events.", call. = FALSE)
		}
	}

	# Appropriate columns for recurrent data
	if (model_type != "trad") {
		if (
			!id %in% names(data) |
			!left %in% names(data) |
			!right %in% names(data) |
			!censor %in% names(data) |
			length(setdiff(event_dates, names(data))) != 0
		) {
			stop("The columns required for recurrent event analysis are not contained with the dataframe.", call. = FALSE)
		}
	}

	# Possible strata
	strata <- paste0("strata_", 0:length(event_dates))

	if (model_type == "trad") {
		tbl <-
			data %>%
			select(all_of(c(id, left, right, censor, event_dates))) %>%
			dplyr::rename(
				id = all_of(id),
				left = all_of(left),
				right = all_of(right),
				censor = all_of(censor)
			) %>%
			mutate(strata = strata) %>%
			mutate(events = sum(!is.na(dplyr::c_across(all_of(event_dates)))))
	} else {
		# Make base table
		tbl <-
			data %>%
			select(all_of(c(id, left, right, censor, event_dates))) %>%
			dplyr::rename(
				id = all_of(id),
				left = all_of(left),
				right = all_of(right),
				censor = all_of(censor)
			) %>%
			dplyr::rowwise() %>%
			mutate(strata_0 = left) %>%
			# Arrange by date
			pivot_longer(
				cols = c(strata_0, all_of(event_dates)),
				names_to = "strata",
				values_to = "date"
			) %>%
			arrange(date) %>%
			# Check all ordered events and rename them (removing same day events)
			group_by(id) %>%
			nest(nested = c(strata, date)) %>%
			mutate(nested = map(nested, function(x) {
				# NA values get trimmed as side effect (will be recovered in pivot wider)
				y <- x %>% group_by(date) %>% dplyr::slice(1)
				n <- length(y$strata)
				y$strata <- strata[1:n]
				return(y)
			})) %>%
			unnest(cols = c(nested)) %>%
			arrange(id, date) %>%
			ungroup()
	}

	# Initialize basic columns for survival table output
	tbl$status <- tbl$start <- tbl$stop <- 0

	# Make final table
	tbl <- recur_model_type(tbl, model_type)

	# Return
	tbl
}

#' @description Switch method for model types
#' @keywords internal
#' @noRd
recur_model_type <- function(tbl, model_type) {

	# Modeling switch
	switch(
		model_type,
		# Traditional / simple model
		trad = {
			tbl$status <- tbl$censor
			tbl$stop <- tbl$right - tbl$left
			tbl$date <- tbl$right
			res <- tbl
		},
		# Marginal Model
		marginal = {
			res <-
				tbl %>%
				group_by(id) %>%
				nest() %>%
				mutate(data = map(data, function(x) {
					# Remove missing rows as they have no data
					x <- na.omit(x)

					# Total number of events, excluding censoring events
					x$events <- nrow(x[x$strata != "strata_0",])
					n <- seq(1:max(x$events))

					# Stop time for events
					for (i in n) {
						# Stop time will be event date - left
						x$stop[x$strata == paste0("strata_", i)] <-
							x$date[x$strata == paste0("strata_", i)] -
							x$left[x$strata == paste0("strata_", i - 1)]

						# Status if event occurs
						x$status[x$strata == paste0("strata_", i)] <- 1
					}

					# Set censor date
					x$date[x$strata == "strata_0"] <-
						x$right[x$strata == "strata_0"]

					# Censoring stop time
					x$stop[x$strata == "strata_0"] <-
						x$date[x$strata == "strata_0"] - x$left[x$strata == "strata_0"]

					# Status of censoring events
					x$status[x$strata == "strata_0" & x$censor == 1] <- 1

					# Return
					return(x)

				})) %>%
				unnest(cols = c(data))
		},
		# Conditional A / PWP Total Time Model
		pwptt = {
			res <-
				tbl %>%
				group_by(id) %>%
				nest() %>%
				mutate(data = map(data, function(x) {
					# Remove missing rows as they have no data
					x <- na.omit(x)

					# Total number of events, excluding censoring events
					x$events <- nrow(x[x$strata != "strata_0",])
					n <- seq(1:max(x$events))

					# Time for events
					for (i in n) {
						# Stop time will be event date - left
						x$stop[x$strata == paste0("strata_", i)] <-
							x$date[x$strata == paste0("strata_", i)] -
							x$left[x$strata == paste0("strata_", i - 1)]

						# Status if event occurs
						x$status[x$strata == paste0("strata_", i)] <- 1

						# Start time
						x$start[x$strata == paste0("strata_", i)] <-
							x$stop[x$strata == paste0("strata_", i - 1)]
					}

					# Set censor date
					x$date[x$strata == "strata_0"] <-
						x$right[x$strata == "strata_0"]

					# Censoring stop time
					x$stop[x$strata == "strata_0"] <-
						x$date[x$strata == "strata_0"] - x$left[x$strata == "strata_0"]

					# Status of censoring events
					x$status[x$strata == "strata_0" & x$censor == 1] <- 1

					# Censor start time
					x$start[x$strata == "strata_0" & x$events > 0] <-
						x$stop[x$strata == paste0("strata_", max(n))]

					# Return
					return(x)

				})) %>%
				unnest(cols = c(data))
		},
		# Conditional B / PWP Gap Time Model
		pwpgt = {
			res <-
				tbl %>%
				group_by(id) %>%
				nest() %>%
				mutate(data = map(data, function(x) {
					# Remove missing rows as they have no data
					x <- na.omit(x)

					# Total number of events, excluding censoring events
					x$events <- nrow(x[x$strata != "strata_0",])
					n <- seq(1:max(x$events))

					# Time for events
					for (i in n) {
						# Stop time will be event date - left (to be modified / collapsed)
						x$stop[x$strata == paste0("strata_", i)] <-
							x$date[x$strata == paste0("strata_", i)] -
							x$left[x$strata == paste0("strata_", i - 1)]

						# Status if event occurs
						x$status[x$strata == paste0("strata_", i)] <- 1

						# Start time will be from prior event
						x$start[x$strata == paste0("strata_", i)] <-
							x$stop[x$strata == paste0("strata_", i - 1)]
					}

					# Set censor date
					x$date[x$strata == "strata_0"] <-
						x$right[x$strata == "strata_0"]

					# Censoring stop time
					x$stop[x$strata == "strata_0"] <-
						x$date[x$strata == "strata_0"] - x$left[x$strata == "strata_0"]

					# Status of censoring events
					x$status[x$strata == "strata_0" & x$censor == 1] <- 1

					# Censor start time
					x$start[x$strata == "strata_0" & x$events > 0] <-
						x$stop[x$strata == paste0("strata_", max(n))]

					# Conditional B has the "start times" collapse down to zero
					x$stop <- x$stop - x$start
					x$start <- x$start - x$start

					# Return
					return(x)

				})) %>%
				unnest(cols = c(data))
		},
		# Counting / Anderson & Gill Model
		ag = {
			res <-
				tbl %>%
				group_by(id) %>%
				nest() %>%
				mutate(data = map(data, function(x) {
					# Remove missing rows as they have no data
					x <- na.omit(x)

					# Total number of events, excluding censoring events
					x$events <- nrow(x[x$strata != "strata_0",])
					n <- seq(1:max(x$events))

					# Time for events
					for (i in n) {
						# Stop time will be event date - left
						x$stop[x$strata == paste0("strata_", i)] <-
							x$date[x$strata == paste0("strata_", i)] -
							x$left[x$strata == paste0("strata_", i - 1)]

						# Status if event occurs
						x$status[x$strata == paste0("strata_", i)] <- 1

						# Start time
						x$start[x$strata == paste0("strata_", i)] <-
							x$stop[x$strata == paste0("strata_", i - 1)]
					}

					# Set censor date
					x$date[x$strata == "strata_0"] <-
						x$right[x$strata == "strata_0"]

					# Censoring stop time
					x$stop[x$strata == "strata_0"] <-
						x$date[x$strata == "strata_0"] - x$left[x$strata == "strata_0"]

					# Status of censoring events
					x$status[x$strata == "strata_0" & x$censor == 1] <- 1

					# Censor start time
					x$start[x$strata == "strata_0" & x$events > 0] <-
						x$stop[x$strata == paste0("strata_", max(n))]

					# For counting model, similar to Conditional A, however only 1 stratum
					x$strata <- "strata_0"

					# Return
					return(x)

				})) %>%
				unnest(cols = c(data))

		},
		# Error catch
		stop(
			paste0("The model type `", model_type, "` is not currently supported.")
		)
	)

	# Return clean results
	res %>%
		arrange(id, date) %>%
		select(c(id, status, start, stop, strata, date, events)) %>%
		ungroup() %>%
		mutate(
			stop = as.numeric(stop),
			start = as.numeric(start)
		)

}
