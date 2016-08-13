library(jsonlite)
N0 = 35000
.data = fromJSON('./Acar20160203-minimal_1.00e-08.json') %>>%
    purrr::map_df(~{
      xy = .[[1]]
      tibble(nu2b= xy[1], nu2f= xy[2], loglik=.[[2]])
    }) %>>%
    dplyr::mutate_at(vars(nu2b, nu2f), function(x) x * N0) %>>%
    (?.)

.max_lik_params = fromJSON('./popt-Acar20160203-minimal-loadpexh_1.00e-08.json') %>>%
    (tibble(nu2b= .[1], nu2f= .[2])) %>>%
    dplyr::mutate_at(vars(nu2b, nu2f), function(x) x * N0) %>>%
    (?.)

max_loglik = max(.data$loglik)

.p = .data %>>%
dplyr::filter(loglik > 1.5 * max_loglik) %>>%
ggplot(aes(x=nu2b, y=nu2f))+
geom_raster(aes(fill=loglik))+
scale_fill_gradient(name='log likelihood', low='#f0f0f0f0', high='#000000')+
geom_point(data=.max_lik_params, colour='#ffffff', size=4)+
scale_x_log10()+
scale_y_log10()+
coord_cartesian(expand=FALSE)+
labs(x=expression(italic(N)["2b"]), y=expression(italic(N)['2f']))+
theme_bw(base_family='sans')+
theme(panel.grid.minor=element_blank())+
theme(legend.position=c(1, 1), legend.justification=c(1, 1))
.p
ggsave('loglik_grey.png', .p, width=5, height=5)

.data %>>%
    dplyr::filter(nu2f == max(nu2f)) %>>%
    dplyr::arrange(desc(loglik))

.data %>>%
    dplyr::filter(nu2f == min(nu2f)) %>>%
    dplyr::arrange(desc(loglik))

