ESPCA = function(X, k=2, overlap.group, k.group=2, we=0.5, t = 0.1, niter=2, err=0.0001, Num.init=3, w_f=weight_f) {
  cat("ESPCA function input diagnostics:\n")
  cat("Class of X:", class(X), "\n")
  cat("Type of X:", typeof(X), "\n")
  cat("Dimensions of X:", dim(X), "\n")
  cat("Is X numeric:", is.numeric(X), "\n")
  cat("Number of NA values in X:", sum(is.na(X)), "\n")
  cat("Class of overlap.group:", class(overlap.group), "\n")
  cat("Length of overlap.group:", length(overlap.group), "\n")
  
  n = nrow(X); 
  p = ncol(X);
  U = matrix(0, n, k); 
  D = matrix(0, k, k); 
  V = matrix(0, p, k);
  tX = X;
  weight_f = w_f;
  out = rank1.ESPCA(tX, overlap.group, k.group, we, t, niter, err, Num.init, w_f=weight_f);
  U[,1] = out$u; 
  V[,1] = out$v; 
  D[1,1] = out$d;
  
  if(k < 2) return (list(U=U, D=D, V=V));
  prev_weights = list()
  v_sum = rep(0, 313783) 
  for (i in 1:313783) {
    prev_weights[[i]] = V[,1][i]
  }
  
  for (i in 2:k) {
    v_sum = rep(0, 313783)
    for (j in 1:313783) {
      v_sum[j] = prev_weights[[j]]  
    } 
    weight_n <- abs(as.numeric(as.character(v_sum))) 
    w_n = weight_n
    tX = tX - c(out$d) * out$u %*% t(out$v);
    UU = U %*% t(U);
    weight_n = w_n;
    out = cycleFun2(tX, UU, overlap.group, k.group, we, t, niter, err, Num.init, w_f, w_n=w_n);
    U[,i] = out$u; 
    V[,i] = out$v; 
    D[i,i] = out$d;
    prev_weights = V[,i-1]  
  }
  return (list(U=U, D=D, V=V));
}

rank1.ESPCA = function(X, overlap.group, k.group, we, t, niter=2, err=0.0001, Num.init = 3,w_f=weight_f){
  # ??????????Ϣ
  cat("rank1.ESPCA function input diagnostics:\n")
  cat("Dimensions of X:", dim(X), "\n")
  cat("Is X numeric:", is.numeric(X), "\n")
  weight_f = w_f
  
  n = nrow(X); # n is the number of samples
  p = ncol(X); # p is the number of features
  d.opt = -100
  # set five initial point
  for(ii in 1:Num.init){
    we1 = we
    print("rank")
    set.seed(ii*100)
    v0 = matrix(rnorm(p,0,1),ncol=1);v0 = v0/norm(v0,'E')
    u0 = matrix(rnorm(n,0,1),ncol=1);u0 = u0/norm(u0,'E')
    # Iterative algorithm to solve u and v values
    for (i in 1:niter){
      # ??????????Ϣ
      cat("Iteration", i, ":\n")
      cat("Class of X %*% v0:", class(X %*% v0), "\n")
      cat("Type of X %*% v0:", typeof(X %*% v0), "\n")
      w_f = weight_f
      u = u.project2(X%*%v0)#
      v = overlap.group.penalty.first (t(X)%*%u, overlap.group, k.group, we1,w_f)
      if(we1 > 0){
        we1 = we1 -t
      }else{
        we1 = 0
      }
      # Algorithm termination condition norm(matrix(v),"E")
      if ((sqrt(sum((u - u0)^2)) <= err) & (sqrt(sum((v - v0)^2)) <= err)){break}
      else {
        u0 = u;v0 = v}
    }
    d =t(u)%*%X%*%v
    if(d>d.opt){
      d.opt = d
      u.opt = u
      v.opt = v
    }
  }
  return (list(u=u.opt, v=v.opt, d=d.opt))
}

cycleFun2 = function(X, UU, overlap.group, k.group, we, t, niter, err, Num.init, w_f, w_n=w_n) {
  cat("cycleFun2 function input diagnostics:\n")
  cat("Dimensions of X:", dim(X), "\n")
  cat("Dimensions of UU:", dim(UU), "\n")
  
  n = nrow(X); 
  p = ncol(X);
  d.opt = -100;
  for(ii in 1:Num.init) {
    print("cyc")
    we1 = we
    set.seed(ii * 100)
    v0 = matrix(rnorm(p, 0, 1), ncol = 1); v0 = v0 / norm(v0, 'E')
    u0 = matrix(rnorm(n, 0, 1), ncol = 1); u0 = u0 / norm(u0, 'E')
    
    for(i in 1:niter) {
      cat("Iteration", i, ":\n")
      cat("Class of X %*% v0:", class(X %*% v0), "\n")
      cat("Type of X %*% v0:", typeof(X %*% v0), "\n")
      
      u = (diag(n) - UU)%*%(X%*%v0); 
      u = u.project2(u)
      
      v = overlap.group.penalty.second(t(X) %*% u, overlap.group, k.group, we1, w_f,  w_n=w_n)
      if(we1 > 0) {
        we1 = we1 - t
      } else {
        we1 = 0
      }
      
      if ((sqrt(sum((u - u0)^2)) <= err) & (sqrt(sum((v - v0)^2)) <= err)){
        break
      } else {
        u0 = u; 
        v0 = v
      }
    }
    d = t(u) %*% X %*% v
    if (d > d.opt) {
      d.opt = d
      u.opt = u
      v.opt = v
    }
  }
  return (list(u = u.opt, v = v.opt, d = d.opt))
}

u.project2 = function(z) {  
  cat("u.project2 function input diagnostics:\n")
  cat("Class of z:", class(z), "\n")
  cat("Type of z:", typeof(z), "\n")
  cat("Length of z:", length(z), "\n")
  
  if (any(is.na(z))) {
    stop("u.project2: Input vector z contains NA values.")
  }
  
  u = z
  if (sum(u^2) == 0) {
    return(rep(0, length(u)))  # 
  } else {
    u = u / sqrt(sum(u^2))  # 
    return(u)
  }
}

overlap.group.penalty.first <- function(u, overlap.group, k0 = 1.0, we = 0.5, w_f, 
                                        population_size = 20, generations = 3) {
  p <- length(u)
  k1 <- min(ceiling(k0 * 1.2), floor(p*0.1))
  
  
  initialize_population <- function() {
    
    group_scores <- sapply(overlap.group, function(g) 
      sqrt(mean(w_f[g] * u[g]^2)))
    top_groups <- order(group_scores, decreasing=TRUE)[1:k0]
    active_features <- unique(unlist(overlap.group[top_groups]))
    greedy_sol <- numeric(p)
    greedy_sol[active_features] <- u[active_features]
    greedy_sol <- greedy_sol / sqrt(sum(greedy_sol^2))
    
    pop <- matrix(0, population_size, p)
    pop[1:3, ] <- t(replicate(3, greedy_sol))
    

    feat_probs <- (w_f^2 + 0.1)/sum(w_f^2 + 0.1)
    for(i in 4:population_size) {
      selected <- sample(p, min(3*k1, p/3), prob=feat_probs)
      pop[i, selected] <- rnorm(length(selected)) * sqrt(w_f[selected])
      pop[i, ] <- pop[i, ] / sqrt(sum(pop[i, ]^2))
    }
    return(pop)
  }
  

  compute_fitness <- function(pop) {
    u_norm <- u / sqrt(sum(u^2))
    projections <- (pop %*% u_norm)^2
    nonzero_counts <- sqrt(rowSums(pop != 0))
    projections - 0.2 * nonzero_counts^1.2
  }
  

  evolve_population <- function(pop, fitness) {

    elite <- pop[which.max(fitness), ]
    

    parents_idx <- sapply(1:population_size, function(i) {
      candidates <- sample(population_size, 3)
      candidates[which.max(fitness[candidates])]
    })
    parents <- pop[parents_idx, ]

    children <- vector("list", population_size/2)
    for(i in 1:(population_size/2)) {
      p1 <- parents[2*i-1, ]
      p2 <- parents[2*i, ]
      
      active_groups <- which(sapply(overlap.group, function(g) 
        any(p1[g]!=0) | any(p2[g]!=0)))
      if(length(active_groups) > 0) {
        cross_group <- sample(active_groups, 1)
        features <- overlap.group[[cross_group]]
        temp <- p1[features]
        p1[features] <- p2[features]
        p2[features] <- temp
      }
      children[[i]] <- list(p1, p2)
    }
    new_pop <- do.call(rbind, unlist(children, recursive=FALSE))
    
   
    mutate <- function(vec) {
      mut_idx <- sample(p, 5, prob = 1/(w_f+0.1))  
      vec[mut_idx] <- vec[mut_idx] + rnorm(5, sd=0.01)
      vec <- vec / sqrt(sum(vec^2))
      vec[abs(vec) < 0.0001] <- 0  
      return(vec)
    }
    new_pop <- t(apply(new_pop, 1, mutate))
    
    # keep elite
    new_pop[1, ] <- elite
    return(new_pop)
  }
  
  # 
  population <- initialize_population()
  for(gen in 1:generations) {
    cat("GA generation", gen, "\n")
    fitness <- compute_fitness(population)
    population <- evolve_population(population, fitness)
  }
  

  final_v <- population[which.max(compute_fitness(population)), ]
  

  if(all(final_v == 0)){
    warning("All-zero vector generated, returning random solution")
    final_v <- rnorm(length(final_v))
  }
  
  final_norm <- sqrt(sum(final_v^2))
  if(final_norm < .Machine$double.eps){
    final_v <- rep(0, length(final_v))
  } else {
    final_v <- final_v / final_norm
  }

  non_zero <- final_v[final_v != 0]
  if(length(non_zero) == 0) return(final_v)
  
  cutoff <- quantile(abs(non_zero), 0.6, na.rm=TRUE)
  final_v[abs(final_v) < cutoff] <- 0
  
  return(final_v)
}
overlap.group.penalty.second <- function(u, overlap.group, k0=1.0, we=0.5, w_f, w_n) {
  p <- length(u)
  
  composite_weights <- 0.3*w_f + 0.7*w_n

  base_sol <- overlap.group.penalty.first(u, overlap.group, k0, we, composite_weights)
  

  active_idx <- which(base_sol != 0)
  

  temporal_boost <- 1 + 0.6*(w_n[active_idx]/max(w_n))
  base_sol[active_idx] <- base_sol[active_idx] * temporal_boost

  cutoff <- median(abs(base_sol[active_idx])) 
  base_sol[abs(base_sol) < cutoff] <- 0
  

  if(all(base_sol == 0)) return(base_sol)
  base_sol / sqrt(sum(base_sol^2))
}
