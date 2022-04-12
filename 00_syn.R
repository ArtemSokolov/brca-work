dl <- synExtra::synDownloader("data")
dl("syn26486688", "syn26541985")

dl_auc <- synExtra::synDownloader("data/aucs")
dl_auc("syn27221487", "syn27221488", "syn27221489")

dl_sig1 <- synExtra::synDownloader("data/sigs/lit")
dl_sig1(synExtra::synGlob("syn26529258", "*.txt"))

dl_sig2 <- synExtra::synDownloader("data/sigs/univ")
dl_sig2(synExtra::synGlob("syn26529816", "*.txt"))