#!/bin/bash

# Define directory structure; for scripts
FIG_SRC=scripts/figures
TAB_SRC=scripts/tables

# Define directory structure; for output
MANUSCRIPT=manuscript
FIG_DIR=manuscript/figures
TAB_DIR=manuscript/tables

# Define directory structure; for logs
LOG=log
LOG_FIG=log/figures
LOG_TAB=log/tables

# Define commands
RFIG=R CMD BATCH --vanilla $< $(LOG_FIG)/$(<F).Rout
RTAB=R CMD BATCH --vanilla $< $(LOG_TAB)/$(<F).Rout

# All
all: ## Run all parts of the makefile
all: figures clean

# Directories
dir_manuscript: ## Make manuscript directory tree
dir_manuscript:
	test ! -d $(MANUSCRIPT) && mkdir $(MANUSCRIPT) || exit 0
	test ! -d $(TAB_DIR) && mkdir $(TAB_DIR) || exit 0
	test ! -d $(FIG_DIR) && mkdir $(FIG_DIR) || exit 0
	
dir_logs: ## Make logs directory tree
dir_logs:
	test ! -d $(LOG) && mkdir $(LOG) || exit 0
	test ! -d $(LOG_FIG) && mkdir $(LOG_FIG) || exit 0
	test ! -d $(LOG_TAB) && mkdir $(LOG_TAB) || exit 0

figures: ## Generate the figures
figures: dir_manuscript \
	dir_logs \
	$(FIG_DIR)/sim_competitive.png \
	$(FIG_DIR)/sim_cooperative.png \
	$(FIG_DIR)/fc_stats.png \
	$(FIG_DIR)/two_factor_ranking.png

# figures
$(FIG_DIR)/sim_competitive.png: $(FIG_SRC)/sim_competitive.R
	$(RFIG)
$(FIG_DIR)/sim_cooperative.png: $(FIG_SRC)/sim_cooperative.R
	$(RFIG)
$(FIG_DIR)/fc_stats.png: $(FIG_SRC)/fc_stats.R
	$(RFIG)
$(FIG_DIR)/two_factor_ranking.png: $(FIG_SRC)/two_factor_ranking.R
	$(RFIG)
	
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