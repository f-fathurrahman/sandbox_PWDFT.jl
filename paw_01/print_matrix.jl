function print_matrix(A, Nrows, Ncols)
    for i in 1:Nrows
        for j in 1:Ncols
            @printf("%10.5f", A[i,j])
        end
        println()
    end
    return
end
