library(pipeR)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

ms_stats = readr::read_tsv('./stats-Acar20160203-minimal-1.00e-08.txt.gz') %>>% (?.)

tidy = ms_stats %>>%
    dplyr::select(-matches('^S_|^Kst$')) %>>%
    tidyr::gather(variable, value) %>>% (?.)

.summary = tidy %>>%
    dplyr::group_by(variable) %>>%
    dplyr::summarise(
        mean=mean(value, na.rm=TRUE),
        sd=sd(value, na.rm=TRUE),
        q001=quantile(value, 0.01, na.rm=TRUE),
        q099=quantile(value, 0.99, na.rm=TRUE)) %>>% (?.)
readr::write_tsv(.summary, 'ms_summary.tsv')

.tr = c(D_1= 'D[Florida]',
        D_2= 'D[Chichi]',
        pi_1= 'pi[Florida]',
        pi_2= 'pi[Chichi]',
        Fst= 'F[ST]')

.labeller = function(labels, multi_line=TRUE) {
    labels$variable = .tr[labels$variable]
    label_parsed(labels, multi_line)
}

.p = tidy %>>% #sample_n(20000) %>>%
    ggplot(aes(value))+
    geom_histogram(bins=50)+
    geom_vline(aes(xintercept=q001), dplyr::filter(.summary, grepl('^D_2', variable)))+
    geom_vline(aes(xintercept=q099), dplyr::filter(.summary, grepl('^pi_2|^Fst', variable)))+
    facet_wrap(~ variable, scale='free_x', ncol=3, labeller=.labeller)+
    theme_bw()+
    theme(panel.grid.minor=element_blank())
.p
ggsave('null_distributions.png', .p, width=7, height=7)
