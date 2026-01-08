library(R6)
library(mlr3) # 0.23 version
library(survival)
library(mlr3proba) # 0.7.1 version
library(mlr3extralearners) # for surv.ranger
library(paradox)
library(survAUC) # 1.4.0
library(mlr3misc)
library(checkmate)
library(survdistr) # 0.0.1@d7babd1

# Define new mlr3(proba) measure that calls `survAUC::AUC.uno()`
# and gets ROC-AUC(t0) at a specific time point t0, using as predictor
# the S(t0), as the measure assumes from the doc "that there is a one-to-one
# relationship between the predictor and the expected survival times conditional
# on the predictor. The predictor here is the survival probability at exactly
# that point.
MeasureSurvUnoAUCt = R6::R6Class("MeasureSurvUnoAUCt",
  inherit = mlr3proba::MeasureSurv,
  public = list(
    initialize = function() {
      param_set = paradox::ps(
        time = p_dbl(0),
        # TRUE => constant S(t) interpolation, otherwise linear interpolation
        constant = p_lgl(init = TRUE)
      )

      super$initialize(
        param_set = param_set,
        id = "surv.uno_auc_t",
        range = c(0, 1),
        minimize = FALSE,
        packages = c("survdistr", "survAUC"),
        predict_type = "distr",
        properties = c("requires_task", "requires_train_set"),
        label = "Uno's ROC-AUC(t)",
        man = "mlr3proba::mlr_measures_surv.uno_auc_t"
      )
    }
  ),

  private = list(
    .score = function(prediction, task, train_set, ...) {
      pv = self$param_set$values
      args = list()

      # time to calculate the AUC(t)
      args$times = checkmate::assert_number(pv$time, na.ok = FALSE, null.ok = FALSE)
      # training data Surv outcome
      args$Surv.rsp = task$truth(train_set)
      # test data Surv outcome
      args$Surv.rsp.new = prediction$truth

      # constant interpolation (default) of S(t) at t = time
      surv = survdistr::mat_interp(
        x = prediction$data$distr, # survival matrix [obs x time points]
        eval_times = args$times,
        constant = pv$constant,
        type = "surv",
        check = FALSE
      )[, 1L]
      # IMPORTANT: need to sign-reverse the S, as survAUC expects predictors
      # where higher values correspond to higher risk, i.e. lower survival
      args$lpnew = -surv
      # cor(prediction$lp, -surv, method = "kend")
      # result is exactly the same as using the linear predictor for cox-like models
      # args$lpnew = prediction$lp

      mlr3misc::invoke(survAUC::AUC.uno, .args = args)$auc
    }
  )
)

# Add measure for easy use with `msr()`
mlr_measures = utils::getFromNamespace("mlr_measures", ns = "mlr3")
mlr_measures$add("surv.uno_auc_t", MeasureSurvUnoAUCt)

# some code to check this measure works as it should
if (FALSE) {
  # get all survival tasks in mlr3proba
  keys = as.data.table(mlr_tasks)[task_type == "surv"][["key"]]
  tasks = lapply(keys, function(key) {
    tsk(key)
  })

  # logging
  lgr::get_logger("mlr3")$set_threshold("warn")

  # Progress bars
  library(progressr)
  options(progressr.enable = TRUE)
  handlers(global = TRUE)
  handlers("progress")

  # parallelization
  future::plan("multisession", workers = 15)

  # conduct benchmark (<=2 min on a modern laptop)
  set.seed(42)
  bm_grid = benchmark_grid(
    tasks = tasks,
    learners = list(
      lrn("surv.kaplan"),
      lrn("surv.ranger", num.trees = 30),
      lrn("surv.coxph")
    ),
    resamplings = rsmp("repeated_cv", repeats = 10, folds = 5)
  )
  bm = benchmark(bm_grid)

  # calculate scores at time ~ median of all tasks time points
  measures = c(
    auc_lp = msr("surv.uno_auc", times = 250), # proba, uses lp
    auc_t  = msr("surv.uno_auc_t", time = 250), # this one, uses S(t)
    auc_t_linear = msr("surv.uno_auc_t", id = "auc_t_linear", time = 250, constant = FALSE) # same as previous, with S(t) linearly interpolated
  )
  res = bm$aggregate(measures)

  # scores don't different almost at all for coxph models => implementation correct
  # we get AUC(t) for RSF as well now!
  res
}
