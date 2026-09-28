library(boinet)

beta_to_pi <- function(beta, S0_tstar, round_digits = 4) {
  # Calculate hazard ratios
  hr <- exp(beta)
  
  # Convert to efficacy probabilities
  pi <- S0_tstar^hr
  
  # Round results
  pi <- round(pi, round_digits)
  
  # Return as named vector if beta has names
  if (!is.null(names(beta))) {
    names(pi) <- names(beta)
  } else {
    names(pi) <- paste0("dose_", 1:length(beta))
  }
  
  return(pi)
}

beta1 <- c(0, -1, -2.7, -3, -2, -1, 0)
beta2 <- c(0, -0.5, -1, -1.5, -2, -2.7, -3)
beta3 <- c(0, -1, -2, -2.7, -2.8, -2.9, -3)
beta4 <- c(0, -0.5, -1, -1.5, -2, -2.3, -3)
beta5 <- c(0, -0.5, -1.0, -1.5, -2.4, -2.7, -3)
beta6 <- c(0, -1.5, -3, -1.5, 0)
beta7 <- c(0, -1, -1.5, -2, -3)
beta8 <- c(0, -1, -2, -2.9, -3)
beta9 <- c(0, -1, -1.5, -1.8, -3)
beta10 <- c(0, -1, -2, -2.7, -3)

pi_1 <- beta_to_pi(beta1, S0_tstar = 0.20)
pi_2 <- beta_to_pi(beta2, S0_tstar = 0.20)
pi_3 <- beta_to_pi(beta3, S0_tstar = 0.20)
pi_4 <- beta_to_pi(beta4, S0_tstar = 0.20)
pi_5 <- beta_to_pi(beta5, S0_tstar = 0.20)
pi_6 <- beta_to_pi(beta6, S0_tstar = 0.20)
pi_7 <- beta_to_pi(beta7, S0_tstar = 0.20)
pi_8 <- beta_to_pi(beta8, S0_tstar = 0.20)
pi_9 <- beta_to_pi(beta9, S0_tstar = 0.20)
pi_10 <- beta_to_pi(beta10, S0_tstar = 0.20)

tox_1 <- c(0.01,0.03,0.05,0.07,0.10,0.12,0.15)
tox_2 <- c(0.01,0.03,0.05,0.07,0.10,0.12,0.15)
tox_3 <- c(0.01,0.03,0.05,0.07,0.10,0.12,0.15)
tox_4 <- c(0.01,0.03,0.05,0.07,0.10,0.15,0.35)
tox_5 <- c(0.01,0.03,0.05,0.07,0.10,0.15,0.35)
tox_6 <- c(0.01,0.03,0.05,0.12,0.15)
tox_7 <- c(0.01,0.03,0.05,0.12,0.15)
tox_8 <- c(0.01,0.03,0.05,0.12,0.15)
tox_9 <- c(0.01,0.03,0.05,0.12,0.35)
tox_10 <- c(0.01,0.03,0.05,0.12,0.35)





tite.boinet(n.dose=7, start.dose=1, size.cohort=4, n.cohort=20,
            toxprob=tox_1, effprob=pi_1,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 80,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=7, start.dose=1, size.cohort=4, n.cohort=20,
            toxprob=tox_2, effprob=pi_2,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 80,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=7, start.dose=1, size.cohort=4, n.cohort=20,
            toxprob=tox_3, effprob=pi_3,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 80,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=7, start.dose=1, size.cohort=4, n.cohort=20,
            toxprob=tox_4, effprob=pi_4,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 80,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=7, start.dose=1, size.cohort=4, n.cohort=20,
            toxprob=tox_5, effprob=pi_5,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 80,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=5, start.dose=1, size.cohort=4, n.cohort=20,
            toxprob=tox_6, effprob=pi_6,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 80,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=5, start.dose=1, size.cohort=4, n.cohort=20,
            toxprob=tox_7, effprob=pi_7,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 80,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=5, start.dose=1, size.cohort=4, n.cohort=20,
            toxprob=tox_8, effprob=pi_8,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 80,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=5, start.dose=1, size.cohort=4, n.cohort=20,
            toxprob=tox_9, effprob=pi_9,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 80,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=5, start.dose=1, size.cohort=4, n.cohort=20,
            toxprob=tox_10, effprob=pi_10,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 80,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)










tite.boinet(n.dose=7, start.dose=1, size.cohort=3, n.cohort=20,
            toxprob=tox_1, effprob=pi_1,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=7, start.dose=1, size.cohort=3, n.cohort=20,
            toxprob=tox_2, effprob=pi_2,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=7, start.dose=1, size.cohort=3, n.cohort=20,
            toxprob=tox_3, effprob=pi_3,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=7, start.dose=1, size.cohort=3, n.cohort=20,
            toxprob=tox_4, effprob=pi_4,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=7, start.dose=1, size.cohort=3, n.cohort=20,
            toxprob=tox_5, effprob=pi_5,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=5, start.dose=1, size.cohort=3, n.cohort=20,
            toxprob=tox_6, effprob=pi_6,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=5, start.dose=1, size.cohort=3, n.cohort=20,
            toxprob=tox_7, effprob=pi_7,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=5, start.dose=1, size.cohort=3, n.cohort=20,
            toxprob=tox_8, effprob=pi_8,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=5, start.dose=1, size.cohort=3, n.cohort=20,
            toxprob=tox_9, effprob=pi_9,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=5, start.dose=1, size.cohort=3, n.cohort=20,
            toxprob=tox_10, effprob=pi_10,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)










tite.boinet(n.dose=7, start.dose=1, size.cohort=4, n.cohort=15,
            toxprob=tox_1, effprob=pi_1,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=7, start.dose=1, size.cohort=4, n.cohort=15,
            toxprob=tox_2, effprob=pi_2,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=7, start.dose=1, size.cohort=4, n.cohort=15,
            toxprob=tox_3, effprob=pi_3,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=7, start.dose=1, size.cohort=4, n.cohort=15,
            toxprob=tox_4, effprob=pi_4,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=7, start.dose=1, size.cohort=4, n.cohort=15,
            toxprob=tox_5, effprob=pi_5,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=5, start.dose=1, size.cohort=4, n.cohort=15,
            toxprob=tox_6, effprob=pi_6,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=5, start.dose=1, size.cohort=4, n.cohort=15,
            toxprob=tox_7, effprob=pi_7,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=5, start.dose=1, size.cohort=4, n.cohort=15,
            toxprob=tox_8, effprob=pi_8,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=5, start.dose=1, size.cohort=4, n.cohort=15,
            toxprob=tox_9, effprob=pi_9,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)


tite.boinet(n.dose=5, start.dose=1, size.cohort=4, n.cohort=15,
            toxprob=tox_10, effprob=pi_10,
            phi = 0.3, phi1 = 0.3*0.1, phi2 = 0.3*1.4,
            delta = 0.9, delta1 = 0.9*0.6,
            alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T = 30, tau.E = 30,
            te.corr = 0.2, gen.event.time = "weibull",
            accrual = 10, gen.enroll.time = "uniform",
            stopping.npts = 60,
            stopping.prob.T = 0.95, stopping.prob.E = 0.99,
            estpt.method = "obs.prob", obd.method = "max.effprob",
            w1 = 0.33, w2 = 1.09,
            plow.ast = 0.3*0.1, pupp.ast = 0.3*1.4,
            qlow.ast = 0.9*0.6/2, qupp.ast = 0.9,
            psi00 = 40, psi11 = 60,
            n.sim = 500, seed.sim = 100)










