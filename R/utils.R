#' Internal parallel lapply wrapper
#'
#' Uses BiocParallel if available and BPPARAM is provided, mclapply on
#' Unix/macOS, or plain lapply on Windows.
#'
#' @param X A vector or list to iterate over.
#' @param FUN The function to apply.
#' @param ... Additional arguments passed to FUN.
#' @param BPPARAM A BiocParallel parameter object, or NULL (default).

robin_lapply <- function(X, FUN, ..., BPPARAM = NULL) {
    if (!is.null(BPPARAM) &&
            requireNamespace("BiocParallel", quietly = TRUE)) {
        BiocParallel::bplapply(X, FUN, ..., BPPARAM = BPPARAM)
    } else if (.Platform$OS.type != "windows") {
        parallel::mclapply(X, FUN, ...)
    } else {
        lapply(X, FUN, ...)
    }
}
