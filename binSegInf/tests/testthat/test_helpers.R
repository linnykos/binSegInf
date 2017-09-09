context("Test helper functions in helpers.R (and elsewhere..).")

test_that("Trimming matrices are done correctly", {


    ## Not used for now, but may be useful for testing in the future..
    alternative_trimrows <- function(mat){
        ind <- apply(mat, 1, function(x) all(is.na(x)))
        return(mat[!ind,,drop=FALSE])
    }

    ## Check row trimming
    mat.list = list(t(matrix(c(1,1,1,1,1,1,NA,NA,NA,NA), nrow=2)),
                    t(matrix(c(1,1,1,1,1,1,1,NA,NA,NA), nrow=2)),
                    t(matrix(c(1,1,NA,1,1,1,NA,NA,NA,NA), nrow=2))
                    )

    trimmed.mat.list = list(t(matrix(c(1,1,1,1,1,1), nrow=2)),
                            t(matrix(c(1,1,1,1,1,1,1,NA), nrow=2)),
                            t(matrix(c(1,1,NA,1,1,1), nrow=2))
                            )

    for(ii in 1:3){
        expect_equal(trim.mat(mat.list[[ii]],type='row'), trimmed.mat.list[[ii]])
    }

    ## Check row+col trimming
    mymat = t(matrix(c(1,NA,1,NA,1,NA,NA,NA,NA,NA), nrow=2))
    mytrimmedmat = (matrix(c(1,1,1), nrow=3))
    expect_equal(trim.mat(mymat, "rowcol"), mytrimmedmat)


})
