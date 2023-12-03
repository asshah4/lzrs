# Test on Yihan's datasets

library(tidyverse)
chf <-
	read_csv("~/Downloads/1. HF_EVENT_final.csv") |>
	janitor::clean_names()

chf <-
	chf |>
	mutate(adhf = case_when(
		chf == "Definite" | chf == "Probable" ~ 1,
		TRUE ~ 0
	)) |>
	mutate(date = as.Date(event_date_4, "%Y/%m/%d")) |>
	rename(id = record_id,
				 lvef = ejection_fract_percent) |>
	select(id, adhf, date, lvef)
