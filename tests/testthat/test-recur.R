# Skip for run time
skip_on_cran()

test_that("recur() has appropriate outputs based on model type", {

	# Setup variables
	library(card)
	data("stress")
	id <- "id"
	first <- "start"
	last <- "stop"
	event_dates <- c(paste0("head_ache_date_", 1:3), paste0("heart_ache_date_", 1:3))
	model_type <- "trad"
	censor <- "death"

	# Traditional
	trad <- recur(stress, model_type = "trad", id, first, last, censor)
	expect_length(trad, 7)
	expect_equal(nrow(trad), nrow(stress))
	expect_equal(trad$strata[1], "strata_0")
	expect_equal(unique(trad$start), 0)

	# Marginal
	marg <- recur(stress, model_type = "marginal", id, first, last, censor, event_dates)
	expect_length(marg, 7)
	expect_gt(nrow(marg), nrow(stress))
	expect_gt(length(unique(marg$strata)), 1)

	# PWP Total Time
	pwptt <- recur(stress, model_type = "pwptt", id, first, last, censor, event_dates)
	expect_gt(nrow(pwptt), nrow(stress))
	expect_gt(length(unique(pwptt$strata)), 1)
	expect_gt(sum(pwptt$start), 0)

	# PWP Gap Time
	pwpgt <- recur(stress, model_type = "pwpgt", id, first, last, censor, event_dates)
	expect_gt(nrow(pwpgt), nrow(stress))
	expect_gt(length(unique(pwpgt$strata)), 1)
	expect_equal(sum(pwpgt$start), 0)

	# Andersen Gill
	ag <- recur(stress, model_type = "ag", id, first, last, censor, event_dates)
	expect_gt(nrow(ag), nrow(stress))
	expect_equal(length(unique(ag$strata)), 1)
	expect_gt(sum(ag$start), 0)

})

