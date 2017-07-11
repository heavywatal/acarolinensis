library(tidyverse)

ms_stats = readr::read_tsv('./stats-Acar20170707-minimal_1.00e-08.tsv.gz') %>% print()

tidy = ms_stats %>%
    dplyr::select(-matches('^S_|^Kst$')) %>%
    tidyr::gather(variable, value) %>%
    tidyr::drop_na() %>%
    print()

.summary = tidy %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
        mean=mean(value, na.rm=TRUE),
        sd=sd(value, na.rm=TRUE),
        q001=quantile(value, 0.01, na.rm=TRUE),
        q005=quantile(value, 0.05, na.rm=TRUE),
        q095=quantile(value, 0.95, na.rm=TRUE),
        q099=quantile(value, 0.99, na.rm=TRUE)) %>% print()
readr::write_tsv(.summary, 'ms_summary.tsv')

.tr = c(D_1= 'D[Florida]',
        D_2= 'D[Ogasawara]',
        pi_1= 'pi[Florida]',
        pi_2= 'pi[Ogasawara]',
        tH_1= 'theta[Florida]',
        tH_2= 'theta[Ogasawara]',
        Fst= 'F[ST]',
        Kst= 'K[ST]')

.labeller = function(labels, multi_line=TRUE) {
    labels$variable = .tr[labels$variable]
    label_parsed(labels, multi_line)
}

.p = tidy %>% #sample_n(20000) %>>%
    dplyr::filter(!str_detect(variable, '_1$')) %>%
    dplyr::filter(!str_detect(variable, '^tH')) %>%
    ggplot(aes(value))+
    geom_histogram(bins=24)+
    geom_point(pch=17, aes(x=q005, y=-2000), dplyr::filter(.summary, grepl('^D_2', variable)))+
    geom_point(pch=17, aes(x=q095, y=-2000), dplyr::filter(.summary, grepl('^pi_2|^Fst', variable)))+
    facet_wrap(~ variable, scales='free_x', ncol=3, labeller=.labeller)+
    theme_bw()+
    theme(axis.title.x=element_blank())+
    theme(panel.grid.minor=element_blank())
.p
ggsave('null_distributions_bin24.png', .p, width=7, height=4)
