#!/bin/bash

# Define directory structure; for scripts
SRC=scripts
FIG_SRC=scripts/figures
TAB_SRC=scripts/tables

# Define directory structure; for output
MANUSCRIPT=manuscript

# Define directory structure; for logs
LOG=log

# Define commands
RFIG=R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
RTAB=R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout

# All
all: ## Run all parts of the makefile
all: figures tables clean

# Directories
dir_manuscript: ## Make manuscript directory tree
dir_manuscript:
	test ! -d $(MANUSCRIPT) && mkdir $(MANUSCRIPT) || exit 0
	
dir_logs: ## Make logs directory tree
dir_logs:
	test ! -d $(LOG) && mkdir $(LOG) || exit 0

figures: ## Generate the figures
figures: dir_manuscript \
	dir_logs \
	$(MANUSCRIPT)/sim_competitive.png \
	$(MANUSCRIPT)/sim_cooperative.png \
	$(MANUSCRIPT)/fc_stats.png \
	$(MANUSCRIPT)/two_factor_ranking.png
	
tables: ## Generate the tables
tables: $(MANUSCRIPT)/tests.tex

# figures
$(MANUSCRIPT)/sim_competitive.png: $(FIG_SRC)/sim_competitive.R
	$(RFIG)
$(MANUSCRIPT)/sim_cooperative.png: $(FIG_SRC)/sim_cooperative.R
	$(RFIG)
$(MANUSCRIPT)/fc_stats.png: $(FIG_SRC)/fc_stats.R $(SRC)/analysis.R
	$(RFIG)
$(MANUSCRIPT)/two_factor_ranking.png: $(FIG_SRC)/two_factor_ranking.R $(SRC)/analysis.R
	$(RFIG)
	
# tables
$(MANUSCRIPT)/tests.tex: $(TAB_SRC)/tests.R $(SRC)/analysis.R
	$(RTAB)
	
# Clean Up
.PHONY: clean
clean: ## Clean up
clean:
	rm -f *.pdf
	rm -f *.RData

# Source: https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
.PHONY: help
help: ## Print the current page
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'
.DEFAULT_GOAL := help