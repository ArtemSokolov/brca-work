dl <- synExtra::synDownloader("data")
dl(
 "syn26486688",     # agents_metadata.csv
 "syn26541985",     # anova.csv
 "syn26560927"      # pmids_for_drugs_tidy.csv
)

dl_auc <- synExtra::synDownloader("data/aucs")
dl_auc("syn27221487", "syn27221488", "syn27221489")
