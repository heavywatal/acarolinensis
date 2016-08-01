library(jsonlite)
N0 = 35000
.data = fromJSON('./Acar20160203-minimal_1.00e-08.json') %>>%
    purrr::map_df(~{
      xy = .[[1]]
      tibble(nu2b= xy[1], nu2f= xy[2], loglik=.[[2]])
    }) %>>% (?.)

.max_lik_params = fromJSON('./popt-Acar20160203-minimal-loadpexh_1.00e-08.json') %>>%
    (tibble(nu2b= .[1], nu2f= .[2])) %>>%
    dplyr::mutate_at(vars(nu2b, nu2f), function(x) x * N0) %>>%
    (?.)

.p = .data %>>%
dplyr::mutate_at(vars(nu2b, nu2f), function(x) x * N0) %>>%
#dplyr::filter(loglik > -1.5e5) %>>%
dplyr::filter(loglik > -3e5) %>>%
# dplyr::filter(loglik > -6e5) %>>%
ggplot(aes(x=nu2b, y=nu2f))+
geom_raster(aes(fill=loglik))+
scale_fill_gradient(low='#f0f0f0f0', high='#404040')+
#scale_fill_distiller(direction=1)+
# geom_contour(aes(z=loglik), colour='#000000', binwidth=4000)+
geom_contour(aes(z=loglik), colour='#000000', binwidth=50000)+
# geom_contour(aes(z=loglik), colour='#000000', size=1, binwidth=100000)+
geom_point(data=.max_lik_params, colour='#000000', size=6)+
scale_x_log10()+
scale_y_log10()+
theme_bw()
.p
ggsave('loglik_contour_grey.png', .p, width=7, height=6)
