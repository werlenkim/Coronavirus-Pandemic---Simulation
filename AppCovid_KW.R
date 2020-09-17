#### R-CODE
#### "CORONAVIRUS PANDEMIC - SIMULATION"
#### AUTHORS: KIM WERLEN & JACOB KOELLA
#### ETH D-USYS 2020

####################################
#### OVERVIEW    

### 1) FUNCTIONS                 
## 1.1) GENERAL FUNCTIONS
## 1.2) TRANSMISSION FUNCTIONS
## 1.3) EPIDEMIOLOGY FUNCTIONS
## 1.4) RATE CALCULATIONS
## 1.5) PLOTS
### 2) INITIALIZATION
## 2.1) GENERAL PARAMETERS
## 2.2) HOMES
## 2.3) AWAY LOCATIONS
## 2.4) HOSPITAL
## 2.5) IDENTIFY PLACES
## 2.6) AGE-ADJUSTED RATES
### 3) MAIN PROGRAM
## 3.1) INITIATION OF POPULATION
## 3.2) INITIATION OF VARIABLES
## 3.3) TIME LOOP
### 4) SHINY APPLICATION
## 4.1) USER INTERFACE
## 4.2) SERVER
## 4.3) CREATING SHINY APP
####################################

rm(list = ls())

# Import libraries
library(shiny)         # required for app
library(matrixStats)   # required for rowProds()
library (fields)       # required for rdist()
library(dplyr)  
library(ggplot2)       # required for ggplot()
library(shinyjs)
library (RANN)         # required for nn2()


### 1) FUNCTIONS ----------------------------------

## 1.1) GENERAL FUNCTIONS ----------------------------------

# going to school or to work

goAway = function(p, closed) {
  if (closed == T) {
    p$x = ifelse(p$xa %in% xpa, p$xa, p$xh)
    p$y = ifelse(p$ya %in% ypa, p$ya, p$yh)
    p$place = ifelse(p$xa %in% xpa, "park", "home")
  } else {
    p$x = p$xa
    p$y = p$ya
    p$place = findPlace(p)
  }
  if (any(p$quarantined)) {
    p$x[p$quarantined] = p$xh[p$quarantined] # set where p$quarantined == TRUE
    p$y[p$quarantined] = p$yh[p$quarantined]
    p$place[p$quarantined] = "home"
  }
  
  return(p)
}

# going home
goHome = function (p) {
  p$x = p$xh
  p$y = p$yh
  if (length(p$place) > 0) {
    p$place = rep("home", length(p$place))
  }
  return (p)
}

# going to the hospital
go2Hospital = function(p) {
  p$x = p$xhos
  p$y = p$yhos
  if (length(p$place) > 0) {
    p$place = rep("hospital", length(p$place))
  }
  return (p)
}

# countdown for Quarantine
countdownQ = function(p) {
  p$timeQ[which(p$timeQ > 0)] = p$timeQ[which(p$timeQ > 0)] - 1 # after a  certain time of isolation, quarantine is over
  # this function is run once per day, so "-1" is substracted
  p$quarantined[which(p$timeQ == 0 &
                        p$quarantined == TRUE)] <- FALSE
  return(p)
}


## 1.2) TRANSMISSION FUNCTION ----------------------------------

# transmission (in person) & contact tracing ----------------------------------
transmission = function(pop, OC, infectMaxA, infectMaxI) {
  S = pop$S
  E = pop$E
  A = pop$A
  I = pop$I
  R = pop$R
  
  if (nrow(I) > 0) {
    nearestNeighbours <- 
      nn2(I[, 2:3],
          # compute distance from infected...
          S[, 2:3],
          k = 1, # we need to know only existence of one neighbour
          searchtype = "radius",
          radius = 0) 
    dI <- nearestNeighbours$nn.dists
    samePlaceI = which(dI == 0) # list of S's indices where distance to I = 0
    s_I  = samePlaceI[runif(length(samePlaceI)) < infectMaxI] # sample percentage that can get infected
    if (length(s_I) > 0) {
      E = rbind(E, S[s_I,])
      S = S[-s_I,]
    }
  } else{
    s_I = c(0)
  } # end of if/else nrow(I)>0
  
  if (nrow(A) > 0) {
    dAA = rdist(A[, 2:3], A[, 2:3]) # A are rows, A are cols
    samePlace = which(dAA == 0, arr.ind = T)
    samePlace[, 1] <- A$ID[samePlace[, 1]]
    samePlace[, 2] <- A$ID[samePlace[, 2]]
    
    # pairwise distance between susceptibles and infecteds
    # also needed for transmission
    if (nrow(S) > 0) {
      nN <- nn2(A[, 2:3],
                S[, 2:3],
                k = 1,
                searchtype = "radius",
                radius = 1)
      dAS <- nN$nn.dists
      samePlaceS_idx = which(dAS == 0)
      
      # transmission
      s_A  = samePlaceS_idx[runif(length(samePlaceS_idx)) < infectMaxA] # sample percentage that can get infected
      if (length(s_A) > 0) {
        E = rbind(E, S[s_A,])
        S = S[-s_A,]
      }
      # contact tracing
      samePlaceS <- matrix(nrow = length(samePlaceS_idx), ncol = 2)
      samePlaceS[, 1] <- A$ID[nN$nn.idx[samePlaceS_idx]]
      samePlaceS[, 2] <- S$ID[samePlaceS_idx]
      if (nrow(samePlace) > 0) {
        samePlace = rbind(samePlaceS, samePlace)
      } else{
        samePlace = samePlaceS
      }
    } else {
      s_A = c(0)
    }
    
    # contact tracing
    if (nrow(E) > 0) {
      nN <- nn2(A[, 2:3],
                E[, 2:3],
                k = 1,
                searchtype = "radius",
                radius = 1)
      dAE <- nN$nn.dists
      samePlaceE_idx = which(dAE == 0)
      samePlaceE = matrix(nrow = length(samePlaceE_idx), ncol = 2)
      samePlaceE[, 1] <- A$ID[nN$nn.idx[samePlaceE_idx]]
      samePlaceE[, 2] <- E$ID[samePlaceE_idx]
      if (nrow(samePlace) > 0) {
        samePlace = rbind(samePlaceE, samePlace)
      } else{
        samePlace = samePlaceE
      }
    }
    # contact tracing
    if (nrow(R) > 0) {
      nN <- nn2(A[, 2:3],
                R[, 2:3],
                k = 1,
                searchtype = "radius",
                radius = 1)
      dAR <- nN$nn.dists
      samePlaceR_idx = which(dAR == 0)
      samePlaceR = matrix(nrow = length(samePlaceR_idx), ncol = 2)
      samePlaceR[, 1] <- A$ID[nN$nn.idx[samePlaceR_idx]]
      samePlaceR[, 2] <- R$ID[samePlaceR_idx]
      if (nrow(samePlace) > 0) {
        samePlace = rbind(samePlaceR, samePlace)
      } else{
        samePlace = samePlaceR
      }
    }
    
    # contact tracing
    spreaders = (unique(samePlace[, 1]))
    if (length(spreaders) > 0) {
      samePlace <- data.frame(x = samePlace[, 1], y = samePlace[, 2])
      contacts <- lapply(1:length(spreaders), function(i) {
        return (c(samePlace[samePlace$x == spreaders[i], 'y']))
      })
    } else{
      contacts <- 0
    }
    
    
  } else{
    contacts <- 0
    spreaders <- 0
    s_A <- c(0)
  } # end of if/else nrow(A) > 0
  
  
  if (typeof(OC) == "list") {
    oldcontacts = OC[[1]]
    oldspreaders = OC[[2]]
    if (typeof(contacts) == "list") {
      contacts <- append(contacts, oldcontacts)
      spreaders <- c(spreaders, oldspreaders)
    } else{
      contacts <- oldcontacts
      spreaders <- oldspreaders
    }
  }
  
  # count successfull transmission events
  sucTrans = c(s_I, s_A)
  
  pop$S = S
  pop$E = E
  pop$A = A
  pop$I = I
  pop$R = R
  return(list(
    pop = pop,
    contacts = contacts,
    spreaders = spreaders,
    sucTrans = sucTrans
  ))
  
}


## 1.3) EPIDEMIOLOGY FUNCTIONS--------------------------

# pass through latent period ----------------------------------
latency = function(pop) {
  E = pop$E
  A = pop$A
  # no differentiation between age classes
  pass = which(sample(
    c(FALSE, TRUE),
    nrow(E),
    replace = T,
    prob = c(1 - latent, latent)
  ))
  if (length(pass) > 0) {
    A = rbind(A, E[pass,])
    E = E[-pass,]
  }
  pop$E = E
  pop$A = A
  return(pop)
}

# develope symptoms & contact tracing ----------------------------------
symptomatic = function(pop) {
  A = pop$A
  I = pop$I
  if (nrow(A) > 0) {
    pass = which(runif(nrow(A), 0, 1) < A$sigma)
    if (length(pass) > 0) {
      # develope symptoms & quarantine
      A$quarantined[pass] <- TRUE
      I = rbind(I, A[pass,])
      A = A[-pass,]
    } else{
      pass = c()
    }
  } else{
    pass = c()
  }
  pop$A = A
  pop$I = I
  return(list(pop = pop, pass = pass))
  
}

# contact tracing ----------------------------------
contactTracing = function(pop, pass, contacts, spreaders, PP) {
  S = pop$S
  E = pop$E
  A = pop$A
  I = pop$I
  R = pop$R
  
  if (typeof(contacts) == "list" & length(pass) > 0) {
    pass.ID = A$ID[pass]
    spreaders.found = which(spreaders %in% pass.ID) # identify position of actual spreaders in vector "spreaders"
    tracing = unlist(contacts[spreaders.found]) # position is the same as in contacts --> find contacts' IDs
    tracing = unique(tracing) # unique(contact's IDs)
    k = runif(length(tracing), 0, 1) < PP # only some participate in contact tracing
    tracing = ifelse(k, tracing, NA)
    
    if (nrow(A) > 0) {
      wa <- which(A$ID %in% tracing) # find contacts by matching IDs
      wa <- wa[!is.na(wa)]
      if (typeof(wa) != "logical") {
        # = logical is the case when wa consists of only NAs --> wa[!is.na(wa)] = logical 0
        A$quarantined[wa] <- TRUE
        A$timeQ[which(A$quarantined == T &
                        A$timeQ == 0)] <-
          isolationTime # those that are newly in quarantine get #isolationTime added
      }
    } else{
      wa = c()
    }
    
    if (nrow(S) > 0) {
      ws <- which(S$ID %in% tracing)
      ws <- ws[!is.na(ws)]
      if (typeof(ws) != "logical") {
        S$quarantined[ws] <- TRUE
        S$timeQ[which(S$quarantined == T &
                        S$timeQ == 0)] <-
          isolationTime # those that are newly in quarantine get #isolationTime added
      }
    } else{
      ws = c()
    }
    
    if (nrow(E) > 0) {
      we <- which(E$ID %in% tracing)
      we <- we[!is.na(we)]
      if (typeof(we) != "logical") {
        E$quarantined[we] <- TRUE
        E$timeQ[which(E$quarantined == T &
                        E$timeQ == 0)] <-
          isolationTime # those that are newly in quarantine get #isolationTime added
      }
    } else{
      we = c()
    }
    
    if (nrow(I) > 0) {
      wi <- which(I$ID %in% tracing)
      wi <- wi[!is.na(wi)]
      if (typeof(wi) != "logical") {
        I$quarantined[wi] <- TRUE
        I$timeQ[which(I$quarantined == T &
                        I$timeQ == 0)] <-
          isolationTime # those that are newly in quarantine get #isolationTime added
      }
    } else{
      wi = c()
    }
    
    if (nrow(R) > 0) {
      wr <- which(R$ID %in% tracing)
      wr <- wr[!is.na(wr)]
      if (typeof(wr) != "logical") {
        R$quarantined[wr] <- TRUE
        R$timeQ[which(R$quarantined == T &
                        R$timeQ == 0)] <-
          isolationTime # those that are newly in quarantine get #isolationTime added
      }
    } else{
      wr = c()
    }
    
    # delete contacts of found spreaders
    contacts[spreaders.found] <- 0
    oldcontacts <- contacts
    oldspreaders <- spreaders
    OC <- list(oldcontacts, oldspreaders)
    quarantineOrders <-
      sum(length(ws), length(wa), length(we), length(wi), length(wr))
    
  } else{
    OC <- 0
    quarantineOrders <- 0
  }
  
  pop$S = S
  pop$E = E
  pop$A = A
  pop$I = I
  pop$R = R
  return(list(pop = pop, OC = OC, qOrder = quarantineOrders))
}

# hospital ----------------------------------
hospitalization = function(pop) {
  I = pop$I
  H = pop$H
  if (nrow(I) > 0) {
    pass = which(runif(nrow(I), 0, 1) < I$hospital)
    if (length(pass) > 0) {
      H = rbind(H, I[pass,])
      I = I[-pass,]
    }
  }
  pop$I = I
  pop$H = H
  return(pop)
}

# death ----------------------------------
death = function(pop) {
  H = pop$H
  D = pop$D
  if (nrow(H) > 0) {
    pass = which(runif(nrow(H), 0, 1) < H$virulence)
    if (length(pass) > 0) {
      D = rbind(D, H[pass,]) 
      H = H[-pass,]
    }
  }
  
  pop$H = H
  pop$D = D
  return(pop)
}



# recovery ----------------------------------
recovery = function(pop) {
  E = pop$E
  A = pop$A
  I = pop$I
  H = pop$H
  R = pop$R
  if (nrow(E) > 0) {
    passE = which(runif(nrow(E), 0, 1) < E$recoveryRateE)
    if (length(passE) > 0) {
      R = rbind(R, E[passE,])
      E = E[-passE,]
    }
  }
  
  if (nrow(A) > 0) {
    passA = which(runif(nrow(A), 0, 1) < A$recoveryRateA)
    if (length(passA) > 0) {
      R = rbind(R, A[passA,])
      A = A[-passA,]
    }
  }
  
  if (nrow(I) > 0) {
    passI = which(runif(nrow(I), 0, 1) < I$recoveryRateI)
    if (length(passI) > 0) {
      R = rbind(R, I[passI,])
      I = I[-passI,]
    }
  }
  
  if (nrow(H) > 0) {
    passH = which(runif(nrow(H), 0, 1) < H$recoveryRateH)
    if (length(passH) > 0) {
      R = rbind(R, H[passH,])
      H = H[-passH,]
    }
  }
  
  pop$R = R
  pop$E = E
  pop$A = A
  pop$H = H
  pop$I = I
  return(pop)
}



## 1.4) RATE CALCULATIONS ----------------------------------

# latency E --> A
getSigma = function(age) {
  # if age ==1 -> 1/2 * (1 + (-1*0.2) =1/2 * (0.8) = 0.4 = 1/2.5 = longer until symptoms develope!
  return(sigma0 * (1 + (age - 2) * 0.2))
}

#recovery E/A/I/H --> R
getRecoveryRate = function(age, comp) {
  # if age ==1 -> 1/2 * (1 - (-1*0.2) =1/2 * (1.2) = 0.6 = 1/1.66 = shorter until recovery!
  return(comp * (1 - (age - 2) * 0.2)) # +
}


## 1.5) PLOTS ----------------------------------

# plot positions of individuals ----------------------------------
drawIndividuals = function(pop, k) {
  plot(
    0,
    type = 'n',
    axes = FALSE,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "",
    ylab = "",
    main = paste("Distribution of People during", k),
    sub = "showing health status"
  )
  box(lwd = 0.8, col = "grey")
  
  lines(
    pop$S$x,
    pop$S$y,
    type = "p",
    col = colors,
    pch = 16,
    cex = 1
  ) 
  lines(
    pop$E$x,
    pop$E$y,
    type = "p",
    col = "yellow2",
    pch = 16,
    cex = 1
  )
  lines(
    pop$A$x,
    pop$A$y,
    type = "p",
    col = "orange",
    pch = 16,
    cex = 1
  )
  lines(
    pop$I$x,
    pop$I$y,
    type = "p",
    col = "red",
    pch = 16,
    cex = 1
  )
  lines(
    pop$H$x,
    pop$H$y,
    type = "p",
    col = "violet",
    pch = 16,
    cex = 1
  )
  lines(
    pop$R$x,
    pop$R$y,
    type = "p",
    col = "green",
    pch = 16,
    cex = 1
  )
  lines(
    pop$D$x,
    pop$D$y,
    type = "p",
    col = "black",
    pch = 16,
    cex = 1
  )
  
  legend(
    "topright",
    legend = c(
      "susceptible",
      "exposed",
      "asymptomatic",
      "symptomatic",
      "hospitalized",
      "recovered",
      "dead"
    ),
    fill = c(
      "dodgerblue",
      "yellow2",
      "orange",
      "red",
      "violet",
      "green",
      "black"
    )
  )
}




# pie charts ----------------------------------
# chart quarantined age prevalence
drawPieQAge = function(pop) {
  slices <-
    c(
      sum(
        length(
          which(pop$S$place == "home" &
                  pop$S$quarantined == T & pop$S$age == 1)
        ),
        length(
          which(pop$E$place == "home" &
                  pop$E$quarantined == T & pop$E$age == 1)
        ),
        length(
          which(pop$A$place == "home" &
                  pop$A$quarantined == T & pop$A$age == 1)
        ),
        length(
          which(pop$I$place == "home" &
                  pop$I$quarantined == T & pop$I$age == 1)
        ),
        length(
          which(pop$R$place == "home" &
                  pop$R$quarantined == T & pop$R$age == 1)
        )
      ),
      
      sum(
        length(
          which(pop$S$place == "home" &
                  pop$S$quarantined == T & pop$S$age == 2)
        ),
        length(
          which(pop$E$place == "home" &
                  pop$E$quarantined == T & pop$E$age == 2)
        ),
        length(
          which(pop$A$place == "home" &
                  pop$A$quarantined == T & pop$A$age == 2)
        ),
        length(
          which(pop$I$place == "home" &
                  pop$I$quarantined == T & pop$I$age == 2)
        ),
        length(
          which(pop$R$place == "home" &
                  pop$R$quarantined == T & pop$R$age == 2)
        )
      ),
      
      sum(
        length(
          which(pop$S$place == "home" &
                  pop$S$quarantined == T & pop$S$age == 3)
        ),
        length(
          which(pop$E$place == "home" &
                  pop$E$quarantined == T & pop$E$age == 3)
        ),
        length(
          which(pop$A$place == "home" &
                  pop$A$quarantined == T & pop$A$age == 3)
        ),
        length(
          which(pop$I$place == "home" &
                  pop$I$quarantined == T & pop$I$age == 3)
        ),
        length(
          which(pop$R$place == "home" &
                  pop$R$quarantined == T & pop$R$age == 3)
        )
      ),
      
      sum(
        length(
          which(pop$S$place == "home" &
                  pop$S$quarantined == F & pop$S$age == 1)
        ),
        length(
          which(pop$E$place == "home" &
                  pop$E$quarantined == F & pop$E$age == 1)
        ),
        length(
          which(pop$A$place == "home" &
                  pop$A$quarantined == F & pop$A$age == 1)
        ),
        length(
          which(pop$I$place == "home" &
                  pop$I$quarantined == F & pop$I$age == 1)
        ),
        length(
          which(pop$R$place == "home" &
                  pop$R$quarantined == F & pop$R$age == 1)
        )
      ),
      
      sum(
        length(
          which(pop$S$place == "home" &
                  pop$S$quarantined == F & pop$S$age == 2)
        ),
        length(
          which(pop$E$place == "home" &
                  pop$E$quarantined == F & pop$E$age == 2)
        ),
        length(
          which(pop$A$place == "home" &
                  pop$A$quarantined == F & pop$A$age == 2)
        ),
        length(
          which(pop$I$place == "home" &
                  pop$I$quarantined == F & pop$I$age == 2)
        ),
        length(
          which(pop$R$place == "home" &
                  pop$R$quarantined == F & pop$R$age == 2)
        )
      ),
      
      sum(
        length(
          which(pop$S$place == "home" &
                  pop$S$quarantined == F & pop$S$age == 3)
        ),
        length(
          which(pop$E$place == "home" &
                  pop$E$quarantined == F & pop$E$age == 3)
        ),
        length(
          which(pop$A$place == "home" &
                  pop$A$quarantined == F & pop$A$age == 3)
        ),
        length(
          which(pop$I$place == "home" &
                  pop$I$quarantined == F & pop$I$age == 3)
        ),
        length(
          which(pop$R$place == "home" &
                  pop$R$quarantined == F & pop$R$age == 3)
        )
      )
    )
  
  lbls <-
    c(
      "Children quarantined",
      "Adults quarantined",
      "Seniors quarantined",
      "Children",
      "Adults",
      "Seniors"
    )
  colos <-
    c("coral",
      "coral3",
      "coral4",
      " darkseagreen1",
      "chartreuse",
      "chartreuse4")
  pct <- round(slices / sum(slices) * 100)
  labls <- paste(pct, "%", sep = "") # ad % to labels
  pie(
    slices,
    labels = labls,
    col = colos,
    radius = 0.5,
    main = paste("Percentage of population in quarantine"),
    sub = "Hospitalized individuals excluded",
    cex.sub = 1.2
  )
  legend(
    "topright",
    title = "Quarantined",
    c("Children", "Adults", "Seniors"),
    cex = 1,
    fill = c("coral", "coral3", "coral4")
  )
  legend(
    "bottomright",
    title = "Not quarantined",
    c("Children", "Adults", "Seniors"),
    cex = 1,
    fill = c(" darkseagreen1", "chartreuse", "chartreuse4")
  )
}

# chart deceased age prevalence
drawPieDead = function (pop) {
  if (nrow(pop$D) > 0) {
    if (length(which(pop$D$age == 1)) > 0) {
      dead_1 <- length(which(pop$D$age == 1))
    } else{
      dead_1 <- 0
    }
    if (length(which(pop$D$age == 2)) > 0) {
      dead_2 <- length(which(pop$D$age == 2))
    } else{
      dead_2 <- 0
    }
    if (length(which(pop$D$age == 3)) > 0) {
      dead_3 <- length(which(pop$D$age == 3))
    } else{
      dead_3 <- 0
    }
    
    slices <- c(dead_1, dead_2, dead_3)
    pct <- round(slices / sum(slices) * 100)
    lbls <- paste(pct, "%", sep = "") # ad % to labels
    pie(
      slices,
      labels = lbls,
      col = c("gray80", "gray50", "black"),
      radius = 0.5,
      main = paste("Age composition of deceased"),
      sub = paste("Number of deceased individuals:", nrow(pop$D)),
      cex.sub = 1.2
    )
    legend(
      "topright",
      c("Children", "Adults", "Seniors"),
      cex = 1,
      fill = c("gray80", "gray50", "black")
    )
  } else{
    legend ("center", legend = "No one has died yet")
  }
}


# plot prevalence ----------------------------------
# calculations for prevalence plot (total population)
calcPrevalence = function(t, pop) {
  totPop = N
  qPop = length(which(pop$S$quarantined == T)) + length(which(pop$E$quarantined ==
                                                                T)) + length(which(pop$A$quarantined == T)) +
    length(which(pop$I$quarantined == T)) + length(which(pop$R$quarantined ==
                                                           T))
  
  pSusceptible[t] <<- nrow(pop$S) / totPop
  pExposed[t] <<- nrow(pop$E) / totPop
  pAsymptomatic[t] <<- nrow(pop$A) / totPop
  pInfected[t] <<- nrow(pop$I) / totPop
  pHospital[t] <<- nrow(pop$H) / totPop
  pRecovered[t] <<- nrow(pop$R) / totPop
  deaths[t] <<- nrow(pop$D) / totPop
  Q[t] <<- qPop / (totPop - nrow(pop$D))
  
  
  return(qPop)
}

# prevalenec plot (total population)
drawPrevalence = function(t) {
  plot(
    1:t,
    pInfected[1:t],
    type = "l",
    col = "red",
    lwd = 2,
    xlim = c(1, t),
    ylim = c(0, 1),
    xlab = "Time (days)",
    ylab = "Prevalence",
    main = " Total Population"
  )
  lines(1:t, pSusceptible[1:t],
        col = "dodgerblue", lwd = 2)
  lines(1:t, pExposed[1:t],
        col = "yellow2", lwd = 2)
  lines(1:t, pAsymptomatic[1:t],
        col = "orange", lwd = 2)
  lines(1:t, pHospital[1:t],
        col = "violet", lwd = 2)
  lines(1:t, pRecovered[1:t],
        col = "green", lwd = 2)
  lines(1:t, deaths[1:t],
        col = "black", lwd = 2)
  lines(1:t, Q[1:t],
        col = "brown", lwd = 2)
  
  legend(
    "topright",
    legend = c(
      "susceptible",
      "exposed",
      "asymptomatic",
      "symptomatic",
      "hospitalized",
      "recovered",
      "dead",
      "quarantined"
    ),
    fill = c(
      "dodgerblue",
      "yellow2",
      "orange",
      "red",
      "violet",
      "green",
      "black",
      "brown"
    )
  )
}

# calculations prevalence plot for different demographics
calcPrevalenceAges = function(t, pop) {
  for (stage in names(pop)) {
    # assign numbers
    assign(paste(stage, "_CH", sep = ""), sum(pop[[stage]]$age == 1))
    assign(paste(stage, "_AD", sep = ""), sum(pop[[stage]]$age == 2))
    assign(paste(stage, "_SE", sep = ""), sum(pop[[stage]]$age == 3))
  }
  
  totPop_CH = S_CH + E_CH + A_CH + I_CH + H_CH + R_CH + D_CH
  totPop_AD = S_AD + E_AD + A_AD + I_AD + H_AD + R_AD + D_AD
  totPop_SE = S_SE + E_SE + A_SE + I_SE + H_SE + R_SE + D_SE
  
  
  pSusceptibleCH[t] <<-   S_CH / totPop_CH
  pExposedCH[t] <<-   E_CH / totPop_CH
  pAsymptomaticCH[t] <<-   A_CH / totPop_CH
  pInfectedCH[t] <<-   I_CH / totPop_CH
  pHospitalCH[t] <<-   H_CH / totPop_CH
  pRecoveredCH[t] <<-   R_CH / totPop_CH
  deathsCH[t] <<-   D_CH / totPop_CH
  
  pSusceptibleAD[t] <<-   S_AD / totPop_AD
  pExposedAD[t] <<-   E_AD / totPop_AD
  pAsymptomaticAD[t] <<- A_AD / totPop_AD
  pInfectedAD[t] <<-   I_AD / totPop_AD
  pHospitalAD[t] <<-   H_AD / totPop_AD
  pRecoveredAD[t] <<-   R_AD / totPop_AD
  deathsAD[t] <<-   D_AD / totPop_AD
  
  pSusceptibleSE[t] <<-   S_SE / totPop_SE
  pExposedSE[t] <<-   E_SE / totPop_SE
  pAsymptomaticSE[t] <<-   A_SE / totPop_SE
  pInfectedSE[t] <<-   I_SE / totPop_SE
  pHospitalSE[t] <<-   H_SE / totPop_SE
  pRecoveredSE[t] <<-   R_SE / totPop_SE
  deathsSE[t] <<-   D_SE / totPop_SE
}

# prevalence plot for children
drawPrevalenceChildren = function(t) {
  plot(
    1:t,
    pInfectedCH[1:t],
    type = "l",
    col = "red",
    lwd = 2,
    xlim = c(1, t),
    ylim = c(0, 1),
    xlab = "Time (days)",
    ylab = "Prevalence",
    main = "Prevalence in Children"
  )
  lines(1:t, pSusceptibleCH[1:t],
        col = "dodgerblue", lwd = 2)
  lines(1:t, pExposedCH[1:t],
        col = "yellow2", lwd = 2)
  lines(1:t, pAsymptomaticCH[1:t],
        col = "orange", lwd = 2)
  lines(1:t, pHospitalCH[1:t],
        col = "violet", lwd = 2)
  lines(1:t, pRecoveredCH[1:t],
        col = "green", lwd = 2)
  lines(1:t, deathsCH[1:t],
        col = "black", lwd = 2)
}

# prevalence plot for adults
drawPrevalenceAdults = function(t) {
  plot(
    1:t,
    pInfectedAD[1:t],
    type = "l",
    col = "red",
    lwd = 2,
    xlim = c(1, t),
    ylim = c(0, 1),
    xlab = "Time (days)",
    ylab = "Prevalence",
    main = "Prevalence in Adults"
  )
  lines(1:t, pSusceptibleAD[1:t],
        col = "dodgerblue", lwd = 2)
  lines(1:t, pExposedAD[1:t],
        col = "yellow2", lwd = 2)
  lines(1:t, pAsymptomaticAD[1:t],
        col = "orange", lwd = 2)
  lines(1:t, pHospitalAD[1:t],
        col = "violet", lwd = 2)
  lines(1:t, pRecoveredAD[1:t],
        col = "green", lwd = 2)
  lines(1:t, deathsAD[1:t],
        col = "black", lwd = 2)
}

# prevalence plot for seniors
drawPrevalenceSeniors = function(t) {
  plot(
    1:t,
    pInfectedSE[1:t],
    type = "l",
    col = "red",
    lwd = 2,
    xlim = c(1, t),
    ylim = c(0, 1),
    xlab = "Time (days)",
    ylab = "Prevalence",
    main = "Prevalence in Seniors"
  )
  lines(1:t, pSusceptibleSE[1:t],
        col = "dodgerblue", lwd = 2)
  lines(1:t, pExposedSE[1:t],
        col = "yellow2", lwd = 2)
  lines(1:t, pAsymptomaticSE[1:t],
        col = "orange", lwd = 2)
  lines(1:t, pHospitalSE[1:t],
        col = "violet", lwd = 2)
  lines(1:t, pRecoveredSE[1:t],
        col = "green", lwd = 2)
  lines(1:t, deathsSE[1:t],
        col = "black", lwd = 2)
}

# scatter plot for places
drawPlaces = function(pop, k) {
  if (nrow(pop$S) > 0) {
    kS <- pop$S[, c("x", "y", "place")]
    kS$comp <- "S"
  } else{
    kS = data.frame(
      x = integer(),
      y = integer(),
      place = character(),
      comp = character()
    )
  }
  if (nrow(pop$E) > 0) {
    kE <- pop$E[, c("x", "y", "place")]
    kE$comp <- "E"
  } else{
    kE = data.frame(
      x = integer(),
      y = integer(),
      place = character(),
      comp = character()
    )
  }
  if (nrow(pop$A) > 0) {
    kA <- pop$A[, c("x", "y", "place")]
    kA$comp <- "A"
  } else{
    kA = data.frame(
      x = integer(),
      y = integer(),
      place = character(),
      comp = character()
    )
  }
  if (nrow(pop$I) > 0) {
    kI <- pop$I[, c("x", "y", "place")]
    kI$comp <- "I"
  } else{
    kI = data.frame(
      x = integer(),
      y = integer(),
      place = character(),
      comp = character()
    )
  }
  if (nrow(pop$H) > 0) {
    kH <- pop$H[, c("x", "y", "place")]
    kH$comp <- "H"
  } else{
    kH = data.frame(
      x = integer(),
      y = integer(),
      place = character(),
      comp = character()
    )
  }
  if (nrow(pop$R) > 0) {
    kR <- pop$R[, c("x", "y", "place")]
    kR$comp <- "R"
  } else{
    kR = data.frame(
      x = integer(),
      y = integer(),
      place = character(),
      comp = character()
    )
  }
  
  collection <- rbind(kS, kE, kA, kI, kH, kR)
  collection$comp <-
    factor(collection$comp, levels = c("S", "E", "A", "I", "H", "R"))
  collection$place <-
    factor(collection$place,
           levels = c("home", "hospital", "park", "school", "work"))
  
  theme_set(theme_gray(base_size = 15))
  
  myColors <-
    c("mediumturquoise",
      "mediumvioletred",
      "seagreen",
      "darkorange",
      "midnightblue")
  names(myColors) <- levels(collection$place)
  colScale <- scale_colour_manual(name = "Place", values = myColors)
  
  p <- ggplot(collection, aes(x, y, color = place)) +
    geom_count() +
    labs(
      subtitle = paste("Distribution of people during the", k),
      caption = paste("Total number of individuals:", nrow(collection)),
      y = "y-coordinate",
      x = "x-coordinate",
      title = "Map of places"
    )
  
  ggPlotPlace <- p + colScale
  
  myColors <-
    c("blue", "yellow2", "orange", "red", "violet", "green")
  names(myColors) <- levels(collection$comp)
  colScale <-
    scale_colour_manual(name = "Health status", values = myColors)
  
  p <- ggplot(collection, aes(x, y, color = comp)) +
    geom_point() +
    labs(
      subtitle = paste("Distribution of people during the", k),
      caption = paste("Total number of individuals:", nrow(collection)),
      y = "y-coordinate",
      x = "x-coordinate",
      title = "Map showing health status"
    )
  ggPlotComp <- p + colScale
  return(list(ggPlotPlace, ggPlotComp))
}


### 2) INITIALIZATION ----------------------------------
## 2.1) GENERAL PARAMETERS ----------------------------------
initPs <- function (maxT) {
  pInfected <<- numeric(length = maxT)
  pSusceptible <<- numeric(length = maxT)
  pExposed <<- numeric(length = maxT)
  pAsymptomatic <<- numeric(length = maxT)
  pInfected <<- numeric(length = maxT)
  pHospital <<- numeric(length = maxT)
  pRecovered <<- numeric(length = maxT)
  deaths <<- numeric(length = maxT)
  
  pInfectedCH <<- numeric(length = maxT)
  pSusceptibleCH <<- numeric(length = maxT)
  pExposedCH <<- numeric(length = maxT)
  pAsymptomaticCH <<- numeric(length = maxT)
  pInfectedCH <<- numeric(length = maxT)
  pHospitalCH <<- numeric(length = maxT)
  pRecoveredCH <<- numeric(length = maxT)
  deathsCH <<- numeric(length = maxT)
  
  pInfectedSE <<- numeric(length = maxT)
  pSusceptibleSE <<- numeric(length = maxT)
  pExposedSE <<- numeric(length = maxT)
  pAsymptomaticSE <<- numeric(length = maxT)
  pInfectedSE <<- numeric(length = maxT)
  pHospitalSE <<- numeric(length = maxT)
  pRecoveredSE <<- numeric(length = maxT)
  deathsSE <<- numeric(length = maxT)
  
  pInfectedAD <<- numeric(length = maxT)
  pSusceptibleAD <<- numeric(length = maxT)
  pExposedAD <<- numeric(length = maxT)
  pAsymptomaticAD <<- numeric(length = maxT)
  pInfectedAD <<- numeric(length = maxT)
  pHospitalAD <<- numeric(length = maxT)
  pRecoveredAD <<- numeric(length = maxT)
  deathsAD <<- numeric(length = maxT)
  Q <<- numeric(length = maxT)
}


isolationTime = 10 # time in quarantine
OC = 0 # old contacts are 0 to begin with
infectMaxI = 0.8 #  p(infection) if distance = 0 for I
infectMaxA = 0.4 #  p(infection) if distance = 0 for A
infectRandom = 0.001 # probability of infection from random encoutners
latent = 1 / 5.4 # passing through latency 
sigma0 = 1 / 3.4 # developing symptoms
recoveryRateE0 = 1 / 5.4 # recovery rate exposed
recoveryRateA0 = 1 / 3.4 # recovery rate asymptomatic
recoveryRateI0 = 1 / 4.4 # recovery rate infected 
CHR = 0.127 # case hospitalization rate
recoveryRateH0 = 1 / 29.5 # recovery rate hospitalized
HFR = 0.42 # hospitalized fatality rate


## 2.2) HOMES ----------------------------------
N = 5500 # approximate size of population; will change later
householdSize = 2.23
nHouseholds = round(N / householdSize)

# members of households ----------------------------------
pSize = c(0.36, 0.33, 0.13, 0.13, 0.06) # probabilities of household sizes
n = sample(1:5, nHouseholds, replace = T, p = pSize)
household = age = xh = yh = c()
for (h in 1:nHouseholds) {
  # first member is either age 2 or 3
  first = sample(2:3, 1, p = c(1, 0.33))
  # elderly in one- or two-person households
  if (first == 3) {
    age = c(age, rep(3, min(n[h], 2)))
    household = c(household, rep(h, min(n[h], 2)))
  } else {
    # if two-person household, second person has age 2
    age = c(age, rep(2, min(n[h], 2)))
    # other household members are children
    if (n[h] > 2) {
      age = c(age, rep(1, n[h] - 2))
    }
    household = c(household, rep(h, n[h]))
  }
}

# True population size ----------------------------------
N = length(age)     

# positions of households ----------------------------------
xHouse = runif(nHouseholds, 0, 1)
yHouse = runif(nHouseholds, 0, 1)

# coordinates for home ----------------------------------
xh = xHouse[household]
yh = yHouse[household]

## 2.3) "AWAY LOCATIONS" (SCHOOLS, WORKPLACES & PARKS) -------------------------------
# positions of schools ----------------------------------
classSize = 20
nKids = round(length(which(age == 1))) # 85% go to school
nSchools <- ceiling(nKids / classSize)
repeat {
  xs = runif(nSchools, 0, 1)
  ys <- runif(nSchools, 0, 1)
  if (sum(xs %in% xh) + sum(ys %in% yh) == 0) {
    # no overlap with homes
    break
  }
}

# positions of work places ----------------------------------
workPlaceSize = 50
nWorkers = round(length(which(age == 2))) # 85% go to work
nTeachers = nSchools  # nSchools = nTeachers since 1 Teacher per school
nWorkPlaces <- round((nWorkers - nTeachers) / workPlaceSize)
repeat {
  xw = runif(nWorkPlaces, 0, 1)
  yw = runif(nWorkPlaces, 0, 1)
  if (length(intersect(xs, xw)) + length(intersect(ys, yw)) == 0 &
      length(intersect(xh, xw)) + length(intersect(yh, yw)) == 0)
    break
}

# position of Parks ----------------------------------
nSeniors = length(which(age == 3))
nOthers = round(0.15 * length(which(age == 1))) + round(0.15 * length(which(age == 2)))
parkSize = 20
nParks = round((nSeniors + nOthers) / parkSize)
repeat {
  xpa = runif(nParks, 0, 1)
  ypa = runif(nParks, 0, 1)
  if ((sum(xpa %in% xh) + sum(ypa %in% yh) == 0) &&
      # no overlap with homes
      (sum(xpa %in% xs) + sum(ypa %in% ys) == 0) &&
      # no overlap with schools
      (sum(xpa %in% xw) + sum(ypa %in% yw) == 0))
    # no overlap with work place
  {
    break
  }
}


# coordinates for school and work and parks ----------------------------------
xa = xh
ya = yh 

park = sample(1:length(xpa), nSeniors, replace = T)
goToPark = runif(nSeniors, 0, 1) < 0.6
xa[which(age == 3)] = ifelse(goToPark, xpa[park], xh)
ya[which(age == 3)] = ifelse(goToPark, ypa[park], yh)

school = sample(1:length(xs), nKids, replace = T)
goToSchool = runif(nKids, 0, 1) < 0.85 
xa[which(age == 1)] = ifelse(goToSchool, xs[school], xpa)
ya[which(age == 1)] = ifelse(goToSchool, ys[school], ypa)

workPlace = sample(1:length(xw), nWorkers, replace = T)
goToWork = runif(nWorkers, 0, 1) < 0.85
preTeacher <- which(goToWork == F)
beTeacher = sample(preTeacher, nTeachers, replace = F)
beTeacherTF <- preTeacher %in% beTeacher
xa[which(age == 2)] = ifelse(goToWork, xw[workPlace], ifelse(beTeacherTF, xs, xpa)) 
ya[which(age == 2)] = ifelse(goToWork, yw[workPlace], ifelse(beTeacherTF, ys, ypa))


## 2.4) HOSPITAL ----------------------------------
# position of hospital ----------------------------------
repeat {
  xhos = runif(1, 0, 1)
  yhos = runif(1, 0, 1)
  if ((sum(xhos %in% xh) + sum(yhos %in% yh) == 0) &&
      # no overlap with homes
      (sum(xhos %in% xs) + sum(yhos %in% ys) == 0) &&
      # no overlap with schools
      (sum(xhos %in% xw) + sum(yhos %in% yw) == 0) &&
      # no overlap with work places
      (sum(xhos %in% xpa) + sum(yhos %in% ypa) == 0)) {
    # no overlap with parks
    break
  }
}

## 2.5) IDENTIFY PLACES ----------------------------------
place = character(length(xh))

findPlace = function(pop) {
  places <- vector(mode = "character", length = nrow(pop))
  if (nrow(pop) > 0) {
    places = ifelse(pop$x %in% xh, "home",
                    ifelse(pop$x %in% xw,
                           "work",
                           ifelse(
                             pop$x %in% xpa,
                             "park",
                             ifelse(
                               pop$x %in% xhos,
                               "hospital",
                               ifelse(pop$x %in% xs,
                                      "school", pop$place)
                             )
                           )))
    
  }
  
  return(places)
}

## 2.6) AGE-ADJUSTED RATES----------------------------------

sigma = getSigma(age)
recoveryRateE = getRecoveryRate(age, recoveryRateE0)
recoveryRateA = getRecoveryRate(age, recoveryRateA0)
recoveryRateI = getRecoveryRate(age, recoveryRateI0)
recoveryRateH = getRecoveryRate(age, recoveryRateH0)
virulence = (recoveryRateH) * HFR / (1 - HFR)
hospital = (recoveryRateI) * CHR / (1 - CHR)



### 3) MAIN PROGRAM ----------------------------------

timeLoop <- function(PP, closed, output, I0, maxTime, TR) {
  
  ## 3.1) INITIATION OF POPULATION ----------------------------------
  initialStatus = c(rep("-", N - I0), rep("+", I0))
  random = sample(1:length(initialStatus), N)
  
  susceptibles = which(initialStatus[random] == "-")
  infecteds    = which(initialStatus[random] == "+")
  initPs(maxTime)
  
  coord = rep(NA, N - I0)
  xho = rep(xhos, N - I0)
  yho = rep(yhos, N - I0)
  S = data.frame(
    ID = susceptibles,
    x = coord,
    y = coord,
    place = character(N - I0),
    xh = xh[susceptibles],
    yh = yh[susceptibles],
    xa = xa[susceptibles],
    ya = ya[susceptibles],
    xhos = xho,
    yhos = yho,
    age = age[susceptibles],
    timeQ = 0,
    quarantined = F,
    sigma = sigma[susceptibles],
    virulence = virulence[susceptibles],
    hospital = hospital[susceptibles],
    recoveryRateE = recoveryRateE[susceptibles],
    recoveryRateA = recoveryRateA[susceptibles],
    recoveryRateI = recoveryRateI[susceptibles],
    recoveryRateH = recoveryRateH[susceptibles]
  )
  rownames(S) <- susceptibles
  coord = rep(NA, I0)
  xho = rep(xhos, I0)
  yho = rep(yhos, I0)
  I = data.frame(
    ID = infecteds,
    x = coord,
    y = coord,
    place = character(I0),
    xh = xh[infecteds],
    yh = yh[infecteds],
    xa = xa[infecteds],
    ya = ya[infecteds],
    xhos = xho,
    yhos = yho,
    age = age[infecteds],
    timeQ = 0,
    quarantined = F,
    sigma = sigma[infecteds],
    virulence = virulence[infecteds],
    hospital = hospital[infecteds],
    recoveryRateE = recoveryRateE[infecteds],
    recoveryRateA = recoveryRateA[infecteds],
    recoveryRateI = recoveryRateI[infecteds],
    recoveryRateH = recoveryRateH[infecteds]
  )
  rownames(I) <- infecteds 
  
  E = data.frame(
    ID = c(),
    x = c(),
    y = c(),
    place = c(),
    xh = c(),
    yh = c(),
    xa = c(),
    ya = c(),
    xhos = c(),
    yhos = c(),
    age = c(),
    timeQ = c(),
    quarantined = c(),
    sigma = c(),
    virulence = c(),
    hospital = c(),
    recoveryRateE = c(),
    recoveryRateA = c(),
    recoveryRateI = c(),
    recoveryRateH = c()
  )
  A = data.frame(
    ID = c(),
    x = c(),
    y = c(),
    place = c(),
    xh = c(),
    yh = c(),
    xa = c(),
    ya = c(),
    xpa = c(),
    yhos = c(),
    age = c(),
    timeQ = c(),
    quarantined = c(),
    sigma = c(),
    virulence = c(),
    hospital = c(),
    recoveryRateE = c(),
    recoveryRateA = c(),
    recoveryRateI = c(),
    recoveryRateH = c()
  )
  H = data.frame(
    ID = c(),
    x = c(),
    y = c(),
    place = c(),
    xh = c(),
    yh = c(),
    xa = c(),
    ya = c(),
    xpa = c(),
    yhos = c(),
    age = c(),
    timeQ = c(),
    quarantined = c(),
    sigma = c(),
    virulence = c(),
    hospital = c(),
    recoveryRateE = c(),
    recoveryRateA = c(),
    recoveryRateI = c(),
    recoveryRateH = c()
  )
  R = data.frame(
    ID = c(),
    x = c(),
    y = c(),
    place = c(),
    xh = c(),
    yh = c(),
    xa = c(),
    ya = c(),
    xpa = c(),
    yhos = c(),
    age = c(),
    timeQ = c(),
    quarantined = c(),
    sigma = c(),
    virulence = c(),
    hospital = c(),
    recoveryRateE = c(),
    recoveryRateA = c(),
    recoveryRateI = c(),
    recoveryRateH = c()
  )
  D = data.frame(
    ID = c(),
    x = c(),
    y = c(),
    place = c(),
    xh = c(),
    yh = c(),
    xa = c(),
    ya = c(),
    xpa = c(),
    yhos = c(),
    age = c(),
    timeQ = c(),
    quarantined = c(),
    sigma = c(),
    virulence = c(),
    hospital = c(),
    recoveryRateE = c(),
    recoveryRateA = c(),
    recoveryRateI = c(),
    recoveryRateH = c()
  )
  
  
  pop = list(
    S = S,
    E = E,
    A = A,
    I = I,
    H = H,
    R = R,
    D = D
  )
  
  ## 3.2) INITIATION OF VARIABLES ----------------------------------
  pop_vec <- list()
  pop_vec[[1]] <- pop
  
  qPopMax <- numeric(N)
  infections <- numeric(2 * N)
  transmissions <- numeric(2 * N)
  transmissions[1] <- 0
  quarantine <- numeric(2 * N)
  quarantine[1] <- 0
  dead <- numeric(N)
  
  infectMaxA <- (infectMaxA * TR)
  infectMaxI <- (infectMaxI * TR)
  
  closed <- closed
  
  ## 3.3) TIME LOOP ----------------------------------
  
  for (t in 1:maxTime) {
    ### 1st Part - DAY
    ## going away
    
    pop = lapply(pop, goAway, closed = closed)
    pop$H = go2Hospital(pop$H)

    ## Transmission and pathogeny
    x = transmission(pop, OC, infectMaxA, infectMaxI)
    pop = x[[1]]
    contacts = x[[2]]
    spreaders = x[[3]]
    sucTrans = x[[4]]
    pop = latency(pop)
    y = symptomatic(pop)
    pop = y[[1]]
    pass = y[[2]]
    z = contactTracing(pop, pass, contacts, spreaders, PP)
    pop = z[[1]]
    OC = z[[2]]
    qOrders = z[[3]]
    pop = hospitalization(pop)
    pop = recovery(pop)
    pop = death(pop)
    
    ## collecting data
    quarantine[2 * t + 1] <- qOrders
    transmissions[2 * t] <- length(sucTrans[which(sucTrans != 0)])
    infections[2 * t] <- length(pass[which(pass != 0)])
    pop_vec[[2 * (t)]] <- pop

    # ----------------------------------------------------
    ### 2nd Part - NIGHT
    ## go home
    pop = lapply(pop, goHome)
    pop$H = go2Hospital(pop$H)

    ## Transmission and pathogeny
    x = transmission(pop, OC, infectMaxA, infectMaxI)
    pop = x[[1]]
    contacts = x[[2]]
    spreaders = x[[3]]
    sucTrans = x[[4]]
    pop = latency(pop)
    y = symptomatic(pop)
    pop = y[[1]]
    pass = y[[2]]
    z = contactTracing(pop, pass, contacts, spreaders, PP)
    pop = z[[1]]
    OC = z[[2]]
    qOrders = z[[3]]
    pop = hospitalization(pop)
    pop = recovery(pop)
    pop = death(pop)
    
    ## collecting data
    quarantine[2 * t + 1] <- qOrders
    transmissions[2 * t + 1] <-
      length(sucTrans[which(sucTrans != 0)])
    infections[2 * t + 1] <- length(pass[which(pass != 0)])
    pop_vec[[2 * (t) + 1]] <- pop
    # ----------------------------------------------------
    
    ### Once per time loop
    ## Quarantine Countdown
    pop = lapply(pop, countdownQ)
    
    ## calculations for plots
    qPopMax[t] = calcPrevalence(t, pop) #use it twice per time step?
    calcPrevalenceAges(t, pop)
    
    ## collecting data
    dead[t] <- nrow(pop$D)
    
  }
  return(
    list(
      pop_vec = pop_vec,
      qPopMax = qPopMax,
      dead = dead,
      infections = infections,
      transmissions = transmissions,
      quarantine = quarantine
    )
  )
}


### 4) SHINY APPLICATION ----------------------------------

####################################
## 4.1) User interface            ##
####################################
ui <- pageWithSidebar(
  # Page header
  headerPanel('Coronavirus Pandemic - Simulation'),
  
  # Input values
  sidebarPanel(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")
    ),
    
    tags$label(h3('Input parameters')),
    sliderInput(
      "PP",
      label = "Percentage of population participating in contact tracing",
      min = 0,
      max = 1,
      step = 0.01,
      value = 0.16
    ),
    hr(),
    h5(strong(
      "Do you want to close schools and workplaces?"
    )),
    checkboxInput("stayhome", label = "Yes", value = F),
    helpText(
      "People will not go to these places anymore, which will reduce transmission opportunities for the virus."
    ),
    hr(),
    sliderInput(
      "TR",
      "Transmission risk (%)",
      min = 0,
      max = 1,
      step = 0.1,
      value = 1
    ),
    helpText(
      "Transmission risk is the probability for a susceptible individual to become infected with the virus after
             having had contact with an infectious person. Transmission risk decreases when people",
      strong("wear masks,"),
      strong("wash their hands regulary and thoroughly"),
      "and practice",
      strong("social distancing.")
    ),
    hr(),
    sliderInput(
      "inf",
      label = "Number of infected people at the beginning",
      min = 1,
      max = round(N / 4) ,
      step = 1,
      value = 20
    ),
    helpText(
      paste(
        "The total number of individuals is",
        N,
        ".",
        "Setting the silder to its maximum value means
                   ~25% of the population is already infected at the start of the pandemic."
      )
    ),
    hr(),
    sliderInput(
      "maxT",
      "Duration of pandemic (days)",
      min = 5,
      max = 50,
      step = 1,
      value = 20
    ),
    helpText(
      "The longer the pandemic, the more time it takes to run the simulation."
    ),
    hr(),
    actionButton("submitbutton", "Submit",
                 class = "btn btn-primary"),
    helpText("Press here to start/update simulation.")
  ),
  
  mainPanel(
    useShinyjs(),
    h2("Information for the user"),
    br(),
    strong("Introduction"),
    p(
      "This interactive simulation explores the effects of the coronavirus pandemic
    on a theoretical population. In the brief summary below, the most important aspects of the
    simulation will be explained. However, if you want to understand the simulation and it’s outputs
    thoroughly, I suggest you take a look at the PDF file",
      a("here.", href = "https://polybox.ethz.ch/index.php/s/7eQa1FydgEAYla1")
    ),
    br(),
    p(
      "R-Code available", a("on GitHub.", href = "")
    ),
    br(),
    strong("Brief summary"),
    p(
      "This is the simulation of a disease outbreak in a fictional city. The population of that city, consisting of children, adults and seniors,
    moves around in the city and visits different places. During the day, schools, work places and parks are visited while at night, people are in their homes.
    Interactions with household members or strangers in public places are opportunities for the virus to be transmitted.
    Depending on  age, the disease progression is more or less severe. Children tend to be less adversly affected than adults,
    while adults are less adversly affected than seniors."
    ),
    p(
      em("The different health stages:"),
      span("Susceptible (S)", style = "color:navy"),
      " = susceptible to the virus;",
      span("Exposed (E)", style = "color:navy"),
      " = was exposed to the virus, not yet infectious;",
      span("Asymptomatic (A)", style = "color:navy"),
      " = infectious but shows no symptoms;",
      span("Symptomatic (I)", style = "color:navy"),
      " = infectious and shows symptoms;",
      span("Hospitalized (H)", style = "color:navy"),
      "= sever symptoms, moved to place Hospital and isolated;",
      span("Recovered (R)", style = "color:navy"),
      " = recovered and immune;",
      span("Dead (D)", style = "color:navy"),
      " = deceased"
    ),
    p(
      em("Quarantine:"),
      "For 10 days, quarantined individuals are not allowed to leave their homes. Symptomatic individuals as well as their contacts who participated in
    contact tracing are quarantined."
    ),
    br(),
    strong("Your turn"),
    p(
      "You can change certain parameters of the simulation and see how the outcome of the pandemic is affected.
    Simply adjust the widgets in the sidebar on the left and press submit to start the simulation.
    "
    ),
    br(),
    
    h2('Status/Output'),
    # Status/Output Text Box
    verbatimTextOutput('contents'),
    tableOutput('tabledata'),
    # Prediction results table,
    # Input for current time
    sliderInput(
      "currentT",
      "Time in pandemic (days)",
      min = 1,
      max = 50,
      step = 1,
      value = 15
    ),
    helpText(
      'Adjusting this slider will allow you to see what the plots look like at a specific time point in the pandemic.
    Adjusting this slider does NOT change the simulation outcome! (silder will appear when simulation is started)'
    ),
    br(),
    # Plots
    plotOutput('prevelancePlot'),
    splitLayout(
      style = "border: white;",
      cellArgs = list(style = "padding: 6px"),
      cellWidths = c("33%", "33%", "33%"),
      plotOutput('prevalenceChildrenPlot'),
      plotOutput('prevalenceAdultsPlot'),
      plotOutput('prevalenceSeniorsPlot')
    ),
    br(),
    splitLayout(
      style = "border: white;",
      cellArgs = list(style = "padding: 6px"),
      cellWidths = c("50%", "50%"),
      plotOutput('cakePlot'),
      plotOutput('cakePlotDead')
    ),
    br(),
    br(),
    h3(
      'The following plots are mainly here to help you understand the environment in which the simulation takes place:'
    ),
    helpText(
      'Try closing schools and work places and see the difference in the plots below.'
    ),
    br(),
    splitLayout(
      style = "border: white;",
      cellArgs = list(style = "padding: 6px"),
      cellWidths = c("50%", "50%"),
      plotOutput('placePlotDay'),
      plotOutput('placePlotNight')
    ),
    br(),
    splitLayout(
      style = "border: white;",
      cellArgs = list(style = "padding: 6px"),
      cellWidths = c("50%", "50%"),
      plotOutput('positionPlotDay'),
      plotOutput('positionPlotNight')
    )
  )
)



####################################
## 4.2) Server                    ##
####################################
server <- function(input, output, session) {

  ############ reactive input listening ###############
  # Adjust input of current time not to exceed max time
  observe({
    maxT <- input$maxT
    newVal <- min(maxT, input$currentT)
    updateSliderInput(
      session,
      "currentT",
      value = newVal,
      min = 1,
      max = maxT
    )
  })
  
  # Hide the time slider when no results are to be found yet
  observe({
    if (input$submitbutton > 0) {
      shinyjs::show(id = "currentT")
    } else {
      shinyjs::hide(id = "currentT")
    }
  })
  
  ############ reactive data collection ###############
  # the simulation result
  # updates upon change of PP, stayhome, inf, maxT, TR
  
  timeLoopResult <- reactive({
    if (!input$submitbutton) {
      return()
    }
    # the following dependencies on pp, etc. are in "isolate" to prevent auto-update
    # on parameter/slider change
    # The use of set.seed with these parameter ensures we get the same random
    # variables for the same parameters.
    isolate(set.seed(sum(
      c(input$PP, input$stayhome, input$inf, input$maxT, input$TR)
    )))
    lastTimeLoopResult <<-
      isolate(timeLoop(
        input$PP,
        input$stayhome,
        output,
        input$inf,
        input$maxT,
        input$TR
      ))
    return(lastTimeLoopResult)
  })
  
  currentPopDay <- reactive({
    res_time_loop <- timeLoopResult()
    pop <- res_time_loop$pop_vec[[input$currentT * 2]]
    return(pop)
  })
  
  currentPopNight <- reactive({
    res_time_loop <- timeLoopResult()
    pop <- res_time_loop$pop_vec[[input$currentT * 2 + 1]]
    return(pop)
  })
  
  ggPlotDay <- reactive({
    ggPlotDay <- drawPlaces(currentPopDay(), "day")
    return(ggPlotDay)
  })
  
  ggPlotNight <- reactive({
    ggPlotNight <- drawPlaces(currentPopNight(), "night")
    return(ggPlotNight)
  })
  
  ############ output assembly ###############
  # Do get output table content
  datasetInput <- reactive({
    res_time_loop <- timeLoopResult()
    qPopMax <- res_time_loop$qPopMax
    infections <- res_time_loop$infections
    dead <- res_time_loop$dead
    transmissions <- res_time_loop$transmissions
    quarantine <- res_time_loop$quarantine
    
    Output <- data.frame(Statistics = c(
      paste(
        "Number of individuals in population =",
        N,
        "(determined by program)"
      ),
      paste(
        "Total number of symptomatic infections:  ",
        sum(infections),
        " = ",
        round(sum(infections) / N, 3) * 100 ,
        "% of the population"
      ),
      paste(
        "Total number of deaths:  ",
        max(dead) ,
        " = ",
        round(max(dead) / N, 3) * 100 ,
        "% of the population"
      ),
      paste("Case Fatality Rate:  ", ifelse(
        sum(infections) != 0, round(max(dead) / sum(infections), 3) * 100, "-"
      ) , "%"),
      paste(
        "Total number of successful transmission events:",
        sum(transmissions),
        " = ",
        round(sum(transmissions) / N, 3) * 100,
        "% of the population was exposed"
      ),
      paste(
        "Total number of people notified by contact tracing app:",
        sum(quarantine),
        " = ",
        round(sum(quarantine) / N, 3) * 100,
        "% of the population"
      ),
      paste (
        "Maximum number of people quarantined simultaneously:",
        max(qPopMax),
        "=",
        round(max(qPopMax) / N, 3) * 100,
        "% of the population"
      )
    ))
    print(Output)
  })
  
  
  ############ output ###############
  ## plots
  # prevalence plots
  output$prevelancePlot <- renderPlot({
    timeLoopResult()
    if (input$submitbutton > 0) {
      drawPrevalence(input$currentT)
    }
  })
  output$prevalenceChildrenPlot <- renderPlot({
    timeLoopResult()
    if (input$submitbutton > 0) {
      drawPrevalenceChildren(input$currentT)
    }
  })
  output$prevalenceAdultsPlot <- renderPlot({
    timeLoopResult()
    if (input$submitbutton > 0) {
      drawPrevalenceAdults(input$currentT)
    }
  })
  output$prevalenceSeniorsPlot <- renderPlot({
    timeLoopResult()
    if (input$submitbutton > 0) {
      drawPrevalenceSeniors(input$currentT)
    }
  })
  
  # cake plots
  output$cakePlot <- renderPlot({
    if (input$submitbutton > 0) {
      drawPieQAge(currentPopNight())
    }
  })
  output$cakePlotDead <- renderPlot({
    if (input$submitbutton > 0) {
      drawPieDead(currentPopNight())
    }
  })
  
  # ggPlots
  output$positionPlotDay <- renderPlot({
    if (input$submitbutton > 0) {
      d <- ggPlotDay()
      d[[2]]
    }
  })
  output$positionPlotNight <- renderPlot({
    if (input$submitbutton > 0) {
      d <- ggPlotNight()
      d[[2]]
    }
  })
  output$placePlotDay <- renderPlot({
    if (input$submitbutton > 0) {
      d <- ggPlotDay()
      d[[1]]
    }
  })
  output$placePlotNight <- renderPlot({
    if (input$submitbutton > 0) {
      d <- ggPlotNight()
      d[[1]]
    }
  })
  
  
  # Status/Output Text Box
  output$contents <- renderPrint({
    if (input$submitbutton > 0) {
      isolate("Simulation complete. Try new parameter settings and start again.")
    } else {
      return("Server is ready for simulation")
    }
  })
  
  # Prediction results table
  output$tabledata <- renderTable({
    if (input$submitbutton > 0) {
      return(datasetInput())
    }
  })
  
  
}

####################################
## 4.3) Create the shiny app      ##
####################################
shinyApp(ui = ui, server = server)