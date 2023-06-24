
# get model weights -------------------------------------------------------

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
