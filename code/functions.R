
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
                sep = '')
        
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
                sep = '')
        
        m1 <- 
          gam(formula = as.formula(form),
              family = poisson,
              offset = log(Effort_sqkm),
              data = data)
        
        models[x] <- 
          AIC(m1)}
      
      form2 <- 
        paste('countInt~', paste0(pull(base_out, output),collapse = '+'))
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

# unscale funciton
unscale <- 
  function(x, data){
    
    scale_info <- 
      scale_factors %>% 
      filter(variable == x)
    
    funct <- 
      function(y){
        ((y * scale_info$scale) + scale_info$center)}
    
    unscaled <- 
      data  %>%  
      mutate_at(x, funct)
    
    return(unscaled)}
  
