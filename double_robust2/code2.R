fit1 <- speff(cd496 ~ age+wtkg+hemo+homo+drugs+karnof+oprior+preanti+
                race+gender+str2+strat+symptom+cd40+cd420+cd80+cd820+offtrt,
              postrandom=c("cd420","cd820","offtrt"), data=ACTG175, trt.id="treat")

fit2 <- speff(cd496 ~ age+wtkg+hemo+homo+drugs+karnof+oprior+preanti+
                race+gender+str2+strat+symptom+cd40+cd420+I(cd420^2)+cd80+cd820+
                I(cd820^2)+cd420:cd820+offtrt, postrandom=c("cd420","I(cd420^2)",
                                                            "cd820","I(cd820^2)","cd420:cd820","offtrt"), data=ACTG175,
              trt.id="treat")

fit3 <- speff(cd496 ~ age+wtkg+hemo+homo+drugs+karnof+oprior+preanti+
                race+gender+str2+strat+symptom+cd40+cd420+cd80+cd820+offtrt,
              postrandom=c("cd420","cd820","offtrt"), data=ACTG175, trt.id="treat",
              optimal="rsq")

fit1 <- speff(cd496 ~ age+wtkg+hemo+homo+drugs+karnof+oprior+preanti+
                race+gender+str2+strat+symptom+cd40+cd420+cd80+cd820+offtrt,
              postrandom=c("cd420","cd820","offtrt"), data=ACTG175, trt.id="treat")
summary(fit1)
