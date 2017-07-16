library(tidyverse)
library(jsonlite)

prefix = 'Acar20160203'
prefix = 'Acar20170707'
N0 = 35000
.exhaustive = fromJSON(paste0(prefix, '-minimal_1.00e-08.json')) %>%
    purrr::map_df(~{
      xy = .[[1]]
      tibble(nu2b= xy[1], nu2f= xy[2], loglik=.[[2]])
    }) %>%
    dplyr::mutate_at(vars(nu2b, nu2f), function(x) x * N0) %>%
    dplyr::arrange(desc(loglik)) %>%
    print()

max_loglik = max(.exhaustive$loglik)

read_param_json = function(infile) {
    fromJSON(infile) %>%
    {tibble(nu2b= .[1], nu2f= .[2])} %>%
    dplyr::mutate_at(vars(nu2b, nu2f), function(x) x * N0)
}

.max_lik_params_opt =
    paste0('popt-', prefix, '-minimal-loadpexh_1.00e-08.json') %>%
    read_param_json() %>% print()

.p = .exhaustive %>%
    dplyr::filter(loglik > 1.2 * max_loglik) %>%
    ggplot(aes(x=nu2b, y=nu2f))+
    geom_raster(aes(fill=loglik))+
    scale_fill_gradient(name='log likelihood', low='#f0f0f0f0', high='#000000')+
    # geom_point(data=top_n(.exhaustive, 1, wt=loglik), colour='#ff9900', size=4)+
    geom_point(data=.max_lik_params_opt, colour='#ffffff', size=4)+
    scale_x_log10()+
    scale_y_log10()+
    coord_cartesian(expand=FALSE)+
    labs(x=expression(italic(N)["2b"]), y=expression(italic(N)['2f']))+
    theme_bw(base_family='sans')+
    theme(panel.grid.minor=element_blank())+
    theme(legend.position=c(1, 1), legend.justification=c(1, 1))
.p
ggsave(paste0('loglik_grey-', prefix, '.png'), .p, width=4, height=4)
