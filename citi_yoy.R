library(dplyr)
library(tidyr)
library(ggplot2)

################################################################################

# Import the daily demand for every start station and bind them into one
# long table

citi_files <- list.files()[grep("citi_", list.files(), fixed = TRUE)]

citi <- readRDS(citi_files[1])
for(i in 2:length(citi_files[1:54])) {
  message(i)
  citi_tmp <- readRDS(citi_files[i])
  citi <- bind_rows(citi, citi_tmp)
}

length(unique(citi$start_station))

################################################################################

# Expand all dates for every start_station

all_dates <- seq(min(citi$date), max(citi$date), by = 1)
citi_full <- expand.grid(id = sort(unique(citi$start_station)), 
                         date = all_dates)
citi_full <- left_join(citi_full, citi, by = c(id = "start_station",
                                               date = "date"))
citi_full <- citi_full %>%
  arrange(id, date) %>%
  rename(rides = n) %>%
  mutate(rides = ifelse(is.na(rides), 0, rides))


head(citi_full)

citi_monthly <- citi_full %>%
  mutate(year = lubridate::year(date),
         month = lubridate::month(date),
         date = as.Date(paste0(substr(date, 1, 7), "-01"))) %>%
  group_by(id, date, year, month) %>%
  summarize(rides = sum(rides)) %>%
  ungroup

citi_yoy <- citi_monthly %>%
  group_by(id) %>%
  mutate(rides_lag12 = lag(rides, 12),
         yoy = (rides - rides_lag12) / rides_lag12,
         yoy = ifelse(is.na(rides_lag12) | (rides_lag12 == 0), NA, yoy),
         yoy_1 = lag(yoy, 1),
         yoy_12 = lag(yoy, 12)) %>%
  ungroup()

citi_yoy %>%
  filter(id == 360) %>%
  ggplot(aes(x = date, y = yoy)) +
  geom_line()

citi_yoy %>%
  filter(id %in% c(360, 384, 472, 491, 3580, 3654)) %>%
  ggplot(aes(x = date, y = yoy, group = as.character(id),
         color = as.character(id))) +
  geom_line()

citi_yoy %>%
  filter(id %in% c(360, 384, 472, 491, 3580, 3654)) %>%
  ggplot(aes(x = yoy_12, y = yoy, group = as.character(id),
             color = as.character(id))) +
  geom_point()

################################################################################

# data for Stan prep

citi_stan <- citi_yoy %>%
  filter(id %in% c(360, 361, 362, 363, 384, 472, 491, 3580, 3654)) %>%
  #filter(id %in% c(360:370, 380:390, 472)) %>%
  filter(!is.na(yoy)) %>%
  filter(!is.na(yoy_1)) %>%
  filter(rides != 0) %>%
  mutate(yoy_12_available = ifelse(is.na(yoy_12), 0, 1)) %>%
  mutate(yoy_12 = ifelse(is.na(yoy_12), 0, yoy_12))

citi_stan <- citi_yoy %>%
  filter(id %in% c(360:500)) %>%
  filter(!is.na(yoy)) %>%
  filter(!is.na(yoy_1)) %>%
  filter(rides != 0) %>%
  mutate(yoy_12_available = ifelse(is.na(yoy_12), 0, 1)) %>%
  mutate(yoy_12 = ifelse(is.na(yoy_12), 0, yoy_12))

N <- dim(citi_stan)[1]
N_stations <- length(unique(citi_stan$id))
station <- citi_stan$id
yoy <- citi_stan$yoy
yoy_1 <- citi_stan$yoy_1
yoy_12 <- citi_stan$yoy_12
yoy_12_available <- citi_stan$yoy_12_available

stan_data <- list(
  N = N,
  N_stations = N_stations,
  station = as.numeric(as.factor(station)),
  yoy = yoy,
  yoy_1 = yoy_1,
  yoy_12 = yoy_12,
  yoy_12_available = yoy_12_available
)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

yoy_12_stan <- stan_model("~/Dropbox/1DataAnalyses/citibike/yoy_12.stan")
yoy_12_fit <- sampling(yoy_12_stan, 
                     data = stan_data,
                     algorithm = "NUTS", seed = 512)
yoy_12_fit <- vb(yoy_12_stan, 
                 data = stan_data,
                 algorithm = "fullrank", seed = 512)

summary(yoy_12_fit)


yoy_data_fit <- as.data.frame(yoy_12_fit)
yoy_hats <- yoy_data_fit %>%
  select(contains("yoy_hat"))
colMeans(yoy_hats)

citi_stan$yoy_hat <- colMeans(yoy_hats)
citi_stan$yoy_q25 <- apply(yoy_hats, 2, quantile, 0.25)
citi_stan$yoy_q75 <- apply(yoy_hats, 2, quantile, 0.75)

citi_stan %>%
  filter(id == 360) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = yoy_q25, ymax = yoy_q75), fill = "blue", alpha = 0.2) +
  geom_point(aes(y = yoy)) +
  geom_line(aes(y = yoy_hat))

citi_stan %>%
  filter(id == 384) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = yoy_q25, ymax = yoy_q75), fill = "blue", alpha = 0.2) +
  geom_point(aes(y = yoy)) +
  geom_line(aes(y = yoy_hat))

citi_stan %>%
  filter(id == 389) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = yoy_q25, ymax = yoy_q75), fill = "blue", alpha = 0.2) +
  geom_point(aes(y = yoy)) +
  geom_line(aes(y = yoy_hat))

citi_stan %>%
  filter(id == 472) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = yoy_q25, ymax = yoy_q75), fill = "blue", alpha = 0.2) +
  geom_point(aes(y = yoy)) +
  geom_line(aes(y = yoy_hat))



################################################################################

yoy_12_ar_stan <- stan_model("~/Dropbox/1DataAnalyses/citibike/yoy_12_ar.stan")
yoy_12_fit <- sampling(yoy_12_ar_stan, 
                       data = stan_data,
                       algorithm = "NUTS", seed = 512)




yoy_data_fit <- as.data.frame(yoy_12_fit)
yoy_hats <- yoy_data_fit %>%
  select(contains("yoy_hat"))
colMeans(yoy_hats)

bs <- yoy_data_fit %>%
  select(contains("b[")) %>%
  summarize_all(mean)

citi_stan$yoy_hat <- colMeans(yoy_hats)
citi_stan$yoy_q25 <- apply(yoy_hats, 2, quantile, 0.25)
citi_stan$yoy_q75 <- apply(yoy_hats, 2, quantile, 0.75)

citi_stan %>%
  filter(id == 360) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = yoy_q25, ymax = yoy_q75), fill = "blue", alpha = 0.2) +
  geom_point(aes(y = yoy)) +
  geom_line(aes(y = yoy_hat))

citi_stan %>%
  filter(id == 384) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = yoy_q25, ymax = yoy_q75), fill = "blue", alpha = 0.2) +
  geom_point(aes(y = yoy)) +
  geom_line(aes(y = yoy_hat))

citi_stan %>%
  filter(id == 386) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = yoy_q25, ymax = yoy_q75), fill = "blue", alpha = 0.2) +
  geom_point(aes(y = yoy)) +
  geom_line(aes(y = yoy_hat))

citi_stan %>%
  filter(id == 389) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = yoy_q25, ymax = yoy_q75), fill = "blue", alpha = 0.2) +
  geom_point(aes(y = yoy)) +
  geom_line(aes(y = yoy_hat))

citi_stan %>%
  filter(id == 472) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = yoy_q25, ymax = yoy_q75), fill = "blue", alpha = 0.2) +
  geom_point(aes(y = yoy)) +
  geom_line(aes(y = yoy_hat))


(0.08 - 0.00232297) / (1 - 0.700718)

(0.0975 - 0.00008765334) / (1 - 0.6926739)




tmp_id <- 21
yoy_orig <- citi_stan %>% filter(id == 472) %>% pull(yoy)
yoy_fut <- rep(NA, 24)

yoy_fut[1] <- mean(yoy_data_fit$a) + mean(yoy_data_fit$`a_station[21]`) +
  mean(yoy_data_fit$`b[21]`) * yoy_orig[length(yoy_orig)]

for(i in 2:length(yoy_fut)) {
  yoy_fut[i] <- mean(yoy_data_fit$a) + mean(yoy_data_fit$`a_station[21]`) +
    mean(yoy_data_fit$`b[21]`) * yoy_fut[i-1]
}

plot(c(yoy_orig, yoy_fut))




tmp_id <- 1
yoy_orig <- citi_stan %>% filter(id == 360) %>% pull(yoy)
yoy_fut <- rep(NA, 24)

yoy_fut[1] <- mean(yoy_data_fit$a) + mean(yoy_data_fit$`a_station[1]`) +
  mean(yoy_data_fit$`b[1]`) * yoy_orig[length(yoy_orig)]

for(i in 2:length(yoy_fut)) {
  yoy_fut[i] <- mean(yoy_data_fit$a) + mean(yoy_data_fit$`a_station[1]`) +
    mean(yoy_data_fit$`b[1]`) * yoy_fut[i-1]
}

plot(c(yoy_orig, yoy_fut))



summary(yoy_12_fit)





###########################################################################


stan_data <- list(
  N = N,
  N_ids = N_stations,
  id = as.numeric(as.factor(station)),
  yoy = yoy,
  yoy_1 = yoy_1
)

dhyoy <- stan_model("~/Dropbox/1DataAnalyses/citibike/damped_hierarchical_yoy.stan")
dhyoy_fit <- sampling(dhyoy, 
                       data = stan_data,
                       algorithm = "NUTS", seed = 512)

yoy_forecast <- function(object, id, h, stan_data) {
  
  yoy <- stan_data$yoy[stan_data$id == id]
  #max_index <- max(which(stan_data$id == id))
  
  #mu_t_train <- as.vector(rstan::extract(object,
  #                                       paste0("mu_t[", max_index, "]"),
  #                                       permuted = FALSE))
  mu <- as.vector(rstan::extract(object, "mu", permuted = FALSE))
  mu_i <- as.vector(rstan::extract(object,
                                   paste0("mu_i[", id, "]"),
                                   permuted = FALSE))
  alpha_i <- as.vector(rstan::extract(object,
                                   paste0("alpha_i[", id, "]"),
                                   permuted = FALSE))
  sigma_i <- as.vector(rstan::extract(object,
                                      paste0("sigma_i[", id, "]"),
                                      permuted = FALSE))
  
  mu_t_test <- matrix(NA, ncol = h, nrow = length(mu))
  yoy_test <- matrix(NA, ncol = h, nrow = length(mu))
  
  mu_t_test[,1] <- (1 - alpha_i) * (mu + mu_i) + alpha_i * yoy[length(yoy)]
  yoy_test[,1] <- rnorm(length(mu), mu_t_test[,1], sigma_i)
  
  for(i in 2:h) {
    mu_t_test[,i] <- (1 - alpha_i) * (mu + mu_i) + alpha_i * yoy_test[,i-1]
    yoy_test[,i] <- rnorm(length(mu), mu_t_test[,i], sigma_i)
  }
  
  aggregated <- data.frame(
    index = 1:(length(yoy) + h),
    yoy = c(yoy, rep(NA, h)),
    yoy_hat = c(rep(NA, length(yoy)), colMeans(yoy_test)),
    yoy_q05 = c(rep(NA, length(yoy)), apply(yoy_test, 2, quantile, 0.05)),
    yoy_q25 = c(rep(NA, length(yoy)), apply(yoy_test, 2, quantile, 0.25)),
    yoy_q50 = c(rep(NA, length(yoy)), apply(yoy_test, 2, quantile, 0.5)),
    yoy_q75 = c(rep(NA, length(yoy)), apply(yoy_test, 2, quantile, 0.75)),
    yoy_q95 = c(rep(NA, length(yoy)), apply(yoy_test, 2, quantile, 0.95)),
    yoy_path_1 = c(rep(NA, length(yoy)), yoy_test[36,]),
    yoy_path_2 = c(rep(NA, length(yoy)), yoy_test[1502,])
  )

  return(
    list(
      yoy = yoy_test,
      mu_t = mu_t_test,
      aggregated = aggregated
    )
  )
}



dhyoy_fc_1 <- yoy_forecast(dhyoy_fit, 1, 24, stan_data)

ggplot(data = dhyoy_fc_1$aggregated) +
  geom_ribbon(aes(x = index, ymin = yoy_q05, ymax = yoy_q95), alpha = 0.2) +
  geom_ribbon(aes(x = index, ymin = yoy_q25, ymax = yoy_q75), alpha = 0.5) +
  geom_line(aes(x = index, y = yoy_hat)) +
  geom_line(aes(x = index, y = yoy)) +
  geom_line(aes(x = index, y = yoy_path_2), color = "orange") +
  geom_line(aes(x = index, y = yoy_path_1), color = "blue")


dhyoy_fc_1 <- yoy_forecast(dhyoy_fit, 7, 24, stan_data)

ggplot(data = dhyoy_fc_1$aggregated) +
  geom_ribbon(aes(x = index, ymin = yoy_q05, ymax = yoy_q95), alpha = 0.2) +
  geom_ribbon(aes(x = index, ymin = yoy_q25, ymax = yoy_q75), alpha = 0.5) +
  geom_line(aes(x = index, y = yoy_hat)) +
  geom_line(aes(x = index, y = yoy)) +
  geom_line(aes(x = index, y = yoy_path_2), color = "orange") +
  geom_line(aes(x = index, y = yoy_path_1), color = "blue")



yoy_data_fit <- as.data.frame(dhyoy_fit)
yoy_hats <- yoy_data_fit %>%
  select(contains("yoy_hat"))
colMeans(yoy_hats)

bs <- yoy_data_fit %>%
  select(contains("b[")) %>%
  summarize_all(mean)

citi_stan$yoy_hat <- colMeans(yoy_hats)
citi_stan$yoy_q25 <- apply(yoy_hats, 2, quantile, 0.25)
citi_stan$yoy_q75 <- apply(yoy_hats, 2, quantile, 0.75)

citi_stan %>%
  filter(id == 360) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = yoy_q25, ymax = yoy_q75), fill = "blue", alpha = 0.2) +
  geom_point(aes(y = yoy)) +
  geom_line(aes(y = yoy_hat))

citi_stan %>%
  filter(id == 384) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = yoy_q25, ymax = yoy_q75), fill = "blue", alpha = 0.2) +
  geom_point(aes(y = yoy)) +
  geom_line(aes(y = yoy_hat))

citi_stan %>%
  filter(id == 386) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = yoy_q25, ymax = yoy_q75), fill = "blue", alpha = 0.2) +
  geom_point(aes(y = yoy)) +
  geom_line(aes(y = yoy_hat))

citi_stan %>%
  filter(id == 389) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = yoy_q25, ymax = yoy_q75), fill = "blue", alpha = 0.2) +
  geom_point(aes(y = yoy)) +
  geom_line(aes(y = yoy_hat))

citi_stan %>%
  filter(id == 472) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = yoy_q25, ymax = yoy_q75), fill = "blue", alpha = 0.2) +
  geom_point(aes(y = yoy)) +
  geom_line(aes(y = yoy_hat))
