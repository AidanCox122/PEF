

# check for strings -------------------------------------------------------
replace_k3_with_k4_stringr <- function(input_string) {
  # Use str_replace_all to replace "k=3" with "k=4" specifically when it follows "phyto", "sst", "salt", or "temp_sd"
  updated_string <- str_replace_all(input_string, "(s\\((phyto|sst|salt|temp_sd),k=)3(\\))", "\\14\\)")
  return(updated_string)
}

# get interannual model weights -------------------------------------------------------

# create a function to calculate AIC weights for different models
getWeightIANN <- 
  # function will accept inputs where base is a vector of established predictors,
  # and test will be a vector of predictors to include in the next iteration of forward selection
  function(base = NULL, test, species = NULL, training = NULL) {
    if(is.null(training)){
      stop('NO TRAINING DATA PROVIDED')
    }
    if(is.null(species)){
      warning('SPECIES CODE NOT PROVIDED, TRAINING DATA MAY CONTAIN MULTIPLE SPECIES')
    }
    # create a repository
    models <- vector(
      mode = 'list',
      length = length(test)) %>% 
      set_names(unique(test))

    data <-
      training %>% 
      filter(Species_code == species)
    
    # now create a reference model with only the base variables
    if(is.null(base)){
      # if base is null, find AIC for each test model
      for(x in names(models)) {
        form <- 
          paste('countInt~',
                's(',
                x,
                ',k=3)',
                sep = '') %>% 
          replace_k3_with_k4_stringr(.)
        
        m1 <- 
          gam(formula = as.formula(form),
              family = poisson,
              offset = log(Effort_sqkm),
              data = data)
        
        models[x] <- 
          AIC(m1)}
      # if base is null, compare to a model with only a random effect of individual
      m2 <- 
        gam(formula = countInt ~ 1,
            family = poisson,
            offset = log(Effort_sqkm),
            data = data)
      # store the results in the models list
      models$null <- 
        AIC(m2)
    } else{
      # otherwise, if base is not null, build a new formula
      
      # separate base if it has multiple components
      # create a repository
    base_out <- tibble(
      output = c()) 
    # create a formula component for each element of base
    for(i in 1:length(base)) {
      output <- tibble(output = paste0('s(', base[i], ',k=3)', collapse = '', sep = ''))
      base_out <- rbind(base_out, output)}
      
    # find AIC for each test model
      for(x in names(models)) {
        form <- 
          paste('countInt~', paste0(pull(base_out, output), collapse = '', sep = '+'),
                's(',
                x,
                ',k=3)',
                sep = '') %>% 
          replace_k3_with_k4_stringr(.)
        
        m1 <- 
          gam(formula = as.formula(form),
              family = poisson,
              offset = log(Effort_sqkm),
              data = data)
        
        models[x] <- 
          AIC(m1)}
      
      form2 <- 
        paste('countInt~', paste0(pull(base_out, output),collapse = '+')) %>% 
        replace_k3_with_k4_stringr(.)
      
      m2 <- 
        gam(formula = as.formula(form2),
            family = poisson,
            offset = log(Effort_sqkm),
            data = data)
      models$null <- 
        AIC(m2)
    }
    
    # calculate AIC weights
    vec_AIC <- unlist(models)
    dAIC <- vec_AIC - min(vec_AIC)
    AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
    return(AICw)
  }


# logistic regression weights (marine mammals) ------------------------------------------

get_logit <- 
  # function will accept inputs where base is a vector of established predictors,
  # and test will be a vector of predictors to include in the next iteration of forward selection
  function(base = NULL, test, species = NULL, training = NULL) {
    if(is.null(training)){
      stop('NO TRAINING DATA PROVIDED')
    }
    if(is.null(species)){
      warning('SPECIES CODE NOT PROVIDED, TRAINING DATA MAY CONTAIN MULTIPLE SPECIES')
    }
    # create a repository
    models <- vector(
      mode = 'list',
      length = length(test)) %>% 
      set_names(unique(test))
    
    data <-
      training %>% 
      filter(Species_code == species)
    
    # find AIC for each test model
    if(is.null(base)){
      # if base is null, find AIC for each test model
      for(x in names(models)) {
        form <- 
          paste('PresAbs~',
                x,
                sep = '')
        
        m1 <- 
          glm(formula = as.formula(form),
              family = 'binomial',
              data = data)
        
        models[x] <- 
          AIC(m1)}
      # if base is null, compare to a model with only a random effect of individual
      m2 <- 
        glm(formula = PresAbs ~ 1,
            family = 'binomial',
            data = data)
      # store the results in the models list
      models$null <- 
        AIC(m2)
    } else{
      # otherwise, if base is not null, build a new formula
      # find AIC for each test model
      for(x in names(models)) {
        form <- paste('PresAbs~', paste0(base,
                            collapse = '',
                            sep = '+'),
              x,
              sep = '')
        
        m1 <- 
          glm(formula = as.formula(form),
              family = 'binomial',
              data = data)
        
        models[x] <- 
          AIC(m1)}
      
      form2 <- 
        paste('PresAbs~', paste(base,collapse = '+'))
      m2 <- 
        glm(formula = as.formula(form2),
            family = 'binomial',
            data = data)
      models$null <- 
        AIC(m2)
    }
    
    # calculate AIC weights
    vec_AIC <- unlist(models)
    dAIC <- vec_AIC - min(vec_AIC)
    AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
    return(AICw)
  }

# GAM weights (seabirds) ------------------------------------------

get_gam <- 
  # function will accept inputs where base is a vector of established predictors,
  # and test will be a vector of predictors to include in the next iteration of forward selection
  function(base = NULL, test, species = NULL, training = NULL) {
    if(is.null(training)){
      stop('NO TRAINING DATA PROVIDED')
    }
    if(is.null(species)){
      warning('SPECIES CODE NOT PROVIDED, TRAINING DATA MAY CONTAIN MULTIPLE SPECIES')
    }
    # create a repository
    models <- vector(
      mode = 'list',
      length = length(test)) %>% 
      set_names(unique(test))
    
    data <-
      training %>% 
      filter(Species_code == species)
    
    # now create a reference model with only the base variables
    if(is.null(base)){
      # if base is null, find AIC for each test model
      for(x in names(models)) {
        form <- 
          paste('Density~',
                's(',
                x,
                ',k=3)',
                sep = '')
        
        m1 <- 
          gam(formula = as.formula(form),
              family = nb,
              data = data)
        
        models[x] <- 
          AIC(m1)}
      # if base is null, compare to a model with only a random effect of individual
      m2 <- 
        gam(formula = Density ~ 1,
            family = nb,
            data = data)
      # store the results in the models list
      models$null <- 
        AIC(m2)
    } else{
      # otherwise, if base is not null, build a new formula
      
      # separate base if it has multiple components
      # create a repository
      base_out <- tibble(
        output = c()) 
      # create a formula component for each element of base
      for(i in 1:length(base)) {
        output <- tibble(output = paste0('s(', base[i], ',k=3)', collapse = '', sep = ''))
        base_out <- rbind(base_out, output)}
      
      # find AIC for each test model
      for(x in names(models)) {
        form <- 
          paste('Density~', paste0(pull(base_out, output), collapse = '', sep = '+'),
                's(',
                x,
                ',k=3)',
                sep = '')
        
        m1 <- 
          gam(formula = as.formula(form),
              family = nb,
              data = data)
        
        models[x] <- 
          AIC(m1)}
      
      form2 <- 
        paste('Density~', paste0(pull(base_out, output),collapse = '+'))
      m2 <- 
        gam(formula = as.formula(form2),
            family = nb,
            data = data)
      models$null <- 
        AIC(m2)
    }
    # calculate AIC weights
    vec_AIC <- unlist(models)
    dAIC <- vec_AIC - min(vec_AIC)
    AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
    return(AICw)
  }


# Unscale Environmenal Variables ------------------------------------------
vars <- c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt', 'dth')

scale_factors <- tibble()
for(x in vars) {
  scale <- 
    daily_mbm_grid %>% 
    pull(x) %>% 
    attr(., 'scaled:scale')
  center <-
    daily_mbm_grid %>% 
    pull(x) %>% 
    attr(., 'scaled:center')
  # store the scaling data in a tibble
  scale_info <- 
    tibble(
      variable = c(x),
      scale = scale,
      center = center)
  # join to repository
  scale_factors <-
    rbind(scale_factors, scale_info)
  #cleanup
  rm(scale_info)
  print(x)}

unscale <- 
  function(x, data, resolution = 'fine'){
    if(resolution != "coarse" & resolution != "fine"){
      stop('Unknown value specified for model resolution. Please specify either fine or coarse')}
    
    if(resolution == 'coarse'){
      scale_info <- 
        coarse_scale_factors %>% 
        filter(variable == x)
      
      funct <- 
        function(y){
          ((y * scale_info$scale) + scale_info$center)}
      
      unscaled <- 
        data  %>%  
        mutate_at(x, funct)
      
      return(unscaled)
    } else{
      print('Model Resolution Specified as: Fine-Scale')
      scale_info <- 
        scale_factors %>% 
        filter(variable == x)
      
      funct <- 
        function(y){
          ((y * scale_info$scale) + scale_info$center)}
      
      unscaled <- 
        data  %>%  
        mutate_at(x, funct)
      
      return(unscaled)}}


# Generalized Additive Mixed Model Weights --------------------------------

get_gamm <- 
  # function will accept inputs where base is a vector of established predictors,
  # and test will be a vector of predictors to include in the next iteration of forward selection
  function(base = NULL, test, random, species = NULL, training = NULL) {
    if(is.null(training)){
      stop('NO TRAINING DATA PROVIDED')}
    
    if(is.null(random)){
      stop('NO RANDOM VARIABLES PROVIDED, USE get_gam() FUNCTION')}
    
    if(any(!random %in% names(training))) {
      warning('Random variable not among training data')}
    if(is.null(species)){
      warning('SPECIES CODE NOT PROVIDED, TRAINING DATA MAY CONTAIN MULTIPLE SPECIES')
    }
    # create a repository
    models <- vector(
      mode = 'list',
      length = length(test)) %>% 
      set_names(unique(test))
    
    data <-
      training %>% 
      filter(Species_code == species)
    
    # create formula components for random variables
    # create a repository
    random_out <- tibble(
      output = c()) 
    # create a formula component for each element of base
    for(i in 1:length(random)) {
      output <- tibble(output = paste0('s(', random[i], ',bs=', print("'re')"), collapse = '', sep = ''))
      random_out <- rbind(random_out, output)}
    
    # now create a reference model with only the base variables
    if(is.null(base)){
      # if base is null, find AIC for each test model
      for(x in names(models)) {
        form <- 
          paste('Density~',
                's(',
                x,
                ',k=4)+',
                paste(pull(random_out, output), collapse = '+'),
                sep = '')
        
        m1 <- 
          gam(formula = as.formula(form),
              family = nb,
              data = data)
        
        models[x] <- 
          AIC(m1)}
      # if base is null, compare to a model with only random effects
      form2 <- paste0(
        'Density~1+', 
        paste(pull(random_out, output), collapse = '+'))
      
      m2 <- 
        gam(formula = as.formula(form2),
            family = nb,
            data = data)
      # store the results in the models list
      models$null <- 
        AIC(m2)
    } else{
      # otherwise, if base is not null, build a new formula
      
      # separate base if it has multiple components
      # create a repository
      base_out <- tibble(
        output = c()) 
      # create a formula component for each element of base
      for(i in 1:length(base)) {
        output <- tibble(output = paste0('s(', base[i], ',k=4)', collapse = '', sep = ''))
        base_out <- rbind(base_out, output)}
      
      # find AIC for each test model
      for(x in names(models)) {
        form <- 
          paste('Density~', paste0(pull(base_out, output), collapse = '', sep = '+'),
                's(',
                x,
                ',k=4)+',
                paste(pull(random_out, output), collapse = '+'),
                sep = '')
        
        m1 <- 
          gam(formula = as.formula(form),
              family = nb,
              data = data)
        
        models[x] <- 
          AIC(m1)}
      
      form2 <- 
        paste0('Density~', paste0(pull(base_out, output),collapse = '', sep = '+'), paste(pull(random_out, output), collapse = '+'))
      m2 <- 
        gam(formula = as.formula(form2),
            family = nb,
            data = data)
      models$null <- 
        AIC(m2)
    }
    # calculate AIC weights
    vec_AIC <- unlist(models)
    dAIC <- vec_AIC - min(vec_AIC)
    AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
    return(AICw)
  }

# Mixed Logistic Regression Weights -----------------------------------------

get_mixed_logit <- 
  # function will accept inputs where base is a vector of established predictors,
  # and test will be a vector of predictors to include in the next iteration of forward selection
  function(base = NULL, test, species = NULL, random = NULL, training = NULL) {
    if(is.null(training)){
      stop('NO TRAINING DATA PROVIDED')
    }
    if(is.null(species)){
      warning('SPECIES CODE NOT PROVIDED, TRAINING DATA MAY CONTAIN MULTIPLE SPECIES')
    }
    # create a repository
    models <- vector(
      mode = 'list',
      length = length(test)) %>% 
      set_names(unique(test))
    
    data <-
      training %>% 
      filter(Species_code == species)
    
    # create formula components for random variables
    # create a repository
    random_out <- tibble(
      output = c()) 
    # create a formula component for each element of base
    for(i in 1:length(random)) {
      output <- tibble(output = paste0('(1|', random[i], ')', collapse = '', sep = ''))
      random_out <- rbind(random_out, output)}
    
    # find AIC for each test model
    if(is.null(base)){
      # if base is null, find AIC for each test model
      for(x in names(models)) {
        form <- 
          paste('PresAbs~',
                x,
                '+',
                paste(pull(random_out, output), collapse = '+'),
                sep = '')
        
        m1 <- 
          glmer(formula = as.formula(form),
              family = 'binomial',
              data = data)
        
        models[x] <- 
          AIC(m1)}
      
      # formula for null model
      form2 <- 
        paste0(
          'PresAbs~1+',
          paste(pull(random_out, output), collapse = '+')
        )
      # if base is null, compare to a model with only a random effect of individual
      m2 <- 
        glmer(formula = as.formula(form2),
            family = 'binomial',
            data = data)
      # store the results in the models list
      models$null <- 
        AIC(m2)
    } else{
      # otherwise, if base is not null, build a new formula
      # find AIC for each test model
      for(x in names(models)) {
        form <- paste('PresAbs~', paste0(base,
                                         collapse = '',
                                         sep = '+'),
                      x,
                      '+',
                      paste(pull(random_out, output), collapse = '+'),
                      sep = '')
        
        m1 <- 
          glmer(formula = as.formula(form),
              family = 'binomial',
              data = data)
        
        models[x] <- 
          AIC(m1)}
      
      form2 <- 
        paste0('PresAbs~',
              paste0(base,
                    '+',
                    paste(pull(random_out, output), collapse = '+'),
                    collapse = '+'))
      m2 <- 
        glmer(formula = as.formula(form2),
            family = 'binomial',
            data = data)
      models$null <- 
        AIC(m2)
    }
    
    # calculate AIC weights
    vec_AIC <- unlist(models)
    dAIC <- vec_AIC - min(vec_AIC)
    AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
    return(AICw)
  }
