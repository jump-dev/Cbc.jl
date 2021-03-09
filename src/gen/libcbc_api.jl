# Julia wrapper for header: Coin_C_defines.h
# Automatically generated using Clang.jl

# Julia wrapper for header: Cbc_C_Interface.h
# Automatically generated using Clang.jl


function Cbc_getVersion()
    ccall((:Cbc_getVersion, libcbcsolver), Cstring, ())
end

function Cbc_newModel()
    ccall((:Cbc_newModel, libcbcsolver), Ptr{Cbc_Model}, ())
end

function Cbc_setProblemName(model, array)
    ccall((:Cbc_setProblemName, libcbcsolver), Cint, (Ptr{Cbc_Model}, Cstring), model, array)
end

function Cbc_addCol(model, name, lb, ub, obj, isInteger, nz, rows, coefs)
    ccall((:Cbc_addCol, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cstring, Cdouble, Cdouble, Cdouble, UInt8, Cint, Ptr{Cint}, Ptr{Cdouble}), model, name, lb, ub, obj, isInteger, nz, rows, coefs)
end

function Cbc_addRow(model, name, nz, cols, coefs, sense, rhs)
    ccall((:Cbc_addRow, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cstring, Cint, Ptr{Cint}, Ptr{Cdouble}, UInt8, Cdouble), model, name, nz, cols, coefs, sense, rhs)
end

function Cbc_addSOS(model, numRows, rowStarts, colIndices, weights, type)
    ccall((:Cbc_addSOS, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint), model, numRows, rowStarts, colIndices, weights, type)
end

function Cbc_loadProblem(model, numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub)
    ccall((:Cbc_loadProblem, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint, Cint, Ptr{CoinBigIndex}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), model, numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub)
end

function Cbc_setColName(model, iColumn, name)
    ccall((:Cbc_setColName, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint, Cstring), model, iColumn, name)
end

function Cbc_setRowName(model, iRow, name)
    ccall((:Cbc_setRowName, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint, Cstring), model, iRow, name)
end

function Cbc_setObjSense(model, sense)
    ccall((:Cbc_setObjSense, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cdouble), model, sense)
end

function Cbc_setRowLower(model, index, value)
    ccall((:Cbc_setRowLower, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint, Cdouble), model, index, value)
end

function Cbc_setRowUpper(model, index, value)
    ccall((:Cbc_setRowUpper, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint, Cdouble), model, index, value)
end

function Cbc_setObjCoeff(model, index, value)
    ccall((:Cbc_setObjCoeff, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint, Cdouble), model, index, value)
end

function Cbc_setColLower(model, index, value)
    ccall((:Cbc_setColLower, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint, Cdouble), model, index, value)
end

function Cbc_setColUpper(model, index, value)
    ccall((:Cbc_setColUpper, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint, Cdouble), model, index, value)
end

function Cbc_setContinuous(model, iColumn)
    ccall((:Cbc_setContinuous, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint), model, iColumn)
end

function Cbc_setInteger(model, iColumn)
    ccall((:Cbc_setInteger, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint), model, iColumn)
end

function Cbc_deleteModel(model)
    ccall((:Cbc_deleteModel, libcbcsolver), Cvoid, (Ptr{Cbc_Model},), model)
end

function Cbc_setMIPStart(model, count, colNames, colValues)
    ccall((:Cbc_setMIPStart, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint, Ptr{Cstring}, Ptr{Cdouble}), model, count, colNames, colValues)
end

function Cbc_setMIPStartI(model, count, colIdxs, colValues)
    ccall((:Cbc_setMIPStartI, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint, Ptr{Cint}, Ptr{Cdouble}), model, count, colIdxs, colValues)
end

function Cbc_clone(model)
    ccall((:Cbc_clone, libcbcsolver), Ptr{Cbc_Model}, (Ptr{Cbc_Model},), model)
end

function Cbc_problemName(model, maxNumberCharacters, array)
    ccall((:Cbc_problemName, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint, Cstring), model, maxNumberCharacters, array)
end

function Cbc_getNumElements(model)
    ccall((:Cbc_getNumElements, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_getNumCols(model)
    ccall((:Cbc_getNumCols, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_getNumIntegers(model)
    ccall((:Cbc_getNumIntegers, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_getNumRows(model)
    ccall((:Cbc_getNumRows, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_getRowName(model, iRow, name, maxLength)
    ccall((:Cbc_getRowName, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint, Cstring, Cint), model, iRow, name, maxLength)
end

function Cbc_getColName(model, iColumn, name, maxLength)
    ccall((:Cbc_getColName, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint, Cstring, Cint), model, iColumn, name, maxLength)
end

function Cbc_getRowNz(model, row)
    ccall((:Cbc_getRowNz, libcbcsolver), Cint, (Ptr{Cbc_Model}, Cint), model, row)
end

function Cbc_getRowIndices(model, row)
    ccall((:Cbc_getRowIndices, libcbcsolver), Ptr{Cint}, (Ptr{Cbc_Model}, Cint), model, row)
end

function Cbc_getRowCoeffs(model, row)
    ccall((:Cbc_getRowCoeffs, libcbcsolver), Ptr{Cdouble}, (Ptr{Cbc_Model}, Cint), model, row)
end

function Cbc_getColNz(model, col)
    ccall((:Cbc_getColNz, libcbcsolver), Cint, (Ptr{Cbc_Model}, Cint), model, col)
end

function Cbc_getColIndices(model, col)
    ccall((:Cbc_getColIndices, libcbcsolver), Ptr{Cint}, (Ptr{Cbc_Model}, Cint), model, col)
end

function Cbc_getColCoeffs(model, col)
    ccall((:Cbc_getColCoeffs, libcbcsolver), Ptr{Cdouble}, (Ptr{Cbc_Model}, Cint), model, col)
end

function Cbc_getRowRHS(model, row)
    ccall((:Cbc_getRowRHS, libcbcsolver), Cdouble, (Ptr{Cbc_Model}, Cint), model, row)
end

function Cbc_getRowSense(model, row)
    ccall((:Cbc_getRowSense, libcbcsolver), UInt8, (Ptr{Cbc_Model}, Cint), model, row)
end

function Cbc_getObjSense(model)
    ccall((:Cbc_getObjSense, libcbcsolver), Cdouble, (Ptr{Cbc_Model},), model)
end

function Cbc_getRowLower(model)
    ccall((:Cbc_getRowLower, libcbcsolver), Ptr{Cdouble}, (Ptr{Cbc_Model},), model)
end

function Cbc_getRowUpper(model)
    ccall((:Cbc_getRowUpper, libcbcsolver), Ptr{Cdouble}, (Ptr{Cbc_Model},), model)
end

function Cbc_getObjCoefficients(model)
    ccall((:Cbc_getObjCoefficients, libcbcsolver), Ptr{Cdouble}, (Ptr{Cbc_Model},), model)
end

function Cbc_getColLower(model)
    ccall((:Cbc_getColLower, libcbcsolver), Ptr{Cdouble}, (Ptr{Cbc_Model},), model)
end

function Cbc_getColUpper(model)
    ccall((:Cbc_getColUpper, libcbcsolver), Ptr{Cdouble}, (Ptr{Cbc_Model},), model)
end

function Cbc_isInteger(model, i)
    ccall((:Cbc_isInteger, libcbcsolver), Cint, (Ptr{Cbc_Model}, Cint), model, i)
end

function Cbc_readMps(model, filename)
    ccall((:Cbc_readMps, libcbcsolver), Cint, (Ptr{Cbc_Model}, Cstring), model, filename)
end

function Cbc_readLp(model, filename)
    ccall((:Cbc_readLp, libcbcsolver), Cint, (Ptr{Cbc_Model}, Cstring), model, filename)
end

function Cbc_writeMps(model, filename)
    ccall((:Cbc_writeMps, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cstring), model, filename)
end

function Cbc_writeLp(model, filename)
    ccall((:Cbc_writeLp, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cstring), model, filename)
end

function Cbc_setInitialSolution(model, sol)
    ccall((:Cbc_setInitialSolution, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Ptr{Cdouble}), model, sol)
end

function Cbc_getVectorStarts(model)
    ccall((:Cbc_getVectorStarts, libcbcsolver), Ptr{CoinBigIndex}, (Ptr{Cbc_Model},), model)
end

function Cbc_getIndices(model)
    ccall((:Cbc_getIndices, libcbcsolver), Ptr{Cint}, (Ptr{Cbc_Model},), model)
end

function Cbc_getElements(model)
    ccall((:Cbc_getElements, libcbcsolver), Ptr{Cdouble}, (Ptr{Cbc_Model},), model)
end

function Cbc_maxNameLength()
    ccall((:Cbc_maxNameLength, libcbcsolver), Cint, ())
end

function Cbc_printModel(model, argPrefix)
    ccall((:Cbc_printModel, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cstring), model, argPrefix)
end

function Cbc_setParameter(model, name, value)
    ccall((:Cbc_setParameter, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cstring, Cstring), model, name, value)
end

function Cbc_getAllowableGap(model)
    ccall((:Cbc_getAllowableGap, libcbcsolver), Cdouble, (Ptr{Cbc_Model},), model)
end

function Cbc_setAllowableGap(model, allowedGap)
    ccall((:Cbc_setAllowableGap, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cdouble), model, allowedGap)
end

function Cbc_getAllowableFractionGap(model)
    ccall((:Cbc_getAllowableFractionGap, libcbcsolver), Cdouble, (Ptr{Cbc_Model},), model)
end

function Cbc_setAllowableFractionGap(model, allowedFracionGap)
    ccall((:Cbc_setAllowableFractionGap, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cdouble), model, allowedFracionGap)
end

function Cbc_getAllowablePercentageGap(model)
    ccall((:Cbc_getAllowablePercentageGap, libcbcsolver), Cdouble, (Ptr{Cbc_Model},), model)
end

function Cbc_setAllowablePercentageGap(model, allowedPercentageGap)
    ccall((:Cbc_setAllowablePercentageGap, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cdouble), model, allowedPercentageGap)
end

function Cbc_getMaximumSeconds(model)
    ccall((:Cbc_getMaximumSeconds, libcbcsolver), Cdouble, (Ptr{Cbc_Model},), model)
end

function Cbc_setMaximumSeconds(model, maxSeconds)
    ccall((:Cbc_setMaximumSeconds, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cdouble), model, maxSeconds)
end

function Cbc_getMaximumNodes(model)
    ccall((:Cbc_getMaximumNodes, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_setMaximumNodes(model, maxNodes)
    ccall((:Cbc_setMaximumNodes, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint), model, maxNodes)
end

function Cbc_getMaximumSolutions(model)
    ccall((:Cbc_getMaximumSolutions, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_setMaximumSolutions(model, maxSolutions)
    ccall((:Cbc_setMaximumSolutions, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint), model, maxSolutions)
end

function Cbc_getLogLevel(model)
    ccall((:Cbc_getLogLevel, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_setLogLevel(model, logLevel)
    ccall((:Cbc_setLogLevel, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cint), model, logLevel)
end

function Cbc_getCutoff(model)
    ccall((:Cbc_getCutoff, libcbcsolver), Cdouble, (Ptr{Cbc_Model},), model)
end

function Cbc_setCutoff(model, cutoff)
    ccall((:Cbc_setCutoff, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, Cdouble), model, cutoff)
end

function Cbc_registerCallBack(model, userCallBack)
    ccall((:Cbc_registerCallBack, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, cbc_callback), model, userCallBack)
end

function Cbc_clearCallBack(model)
    ccall((:Cbc_clearCallBack, libcbcsolver), Cvoid, (Ptr{Cbc_Model},), model)
end

function Cbc_addCutCallback(model, cutcb, name, appData)
    ccall((:Cbc_addCutCallback, libcbcsolver), Cvoid, (Ptr{Cbc_Model}, cbc_cut_callback, Cstring, Ptr{Cvoid}), model, cutcb, name, appData)
end

function Cbc_solve(model)
    ccall((:Cbc_solve, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_getColSolution(model)
    ccall((:Cbc_getColSolution, libcbcsolver), Ptr{Cdouble}, (Ptr{Cbc_Model},), model)
end

function Cbc_getBestPossibleObjValue(model)
    ccall((:Cbc_getBestPossibleObjValue, libcbcsolver), Cdouble, (Ptr{Cbc_Model},), model)
end

function Cbc_bestSolution(model)
    ccall((:Cbc_bestSolution, libcbcsolver), Ptr{Cdouble}, (Ptr{Cbc_Model},), model)
end

function Cbc_numberSavedSolutions(model)
    ccall((:Cbc_numberSavedSolutions, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_savedSolution(model, whichSol)
    ccall((:Cbc_savedSolution, libcbcsolver), Ptr{Cdouble}, (Ptr{Cbc_Model}, Cint), model, whichSol)
end

function Cbc_savedSolutionObj(model, whichSol)
    ccall((:Cbc_savedSolutionObj, libcbcsolver), Cdouble, (Ptr{Cbc_Model}, Cint), model, whichSol)
end

function Cbc_getReducedCost(model)
    ccall((:Cbc_getReducedCost, libcbcsolver), Ptr{Cdouble}, (Ptr{Cbc_Model},), model)
end

function Cbc_isAbandoned(model)
    ccall((:Cbc_isAbandoned, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_isProvenOptimal(model)
    ccall((:Cbc_isProvenOptimal, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_isProvenInfeasible(model)
    ccall((:Cbc_isProvenInfeasible, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_isContinuousUnbounded(model)
    ccall((:Cbc_isContinuousUnbounded, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_getObjValue(model)
    ccall((:Cbc_getObjValue, libcbcsolver), Cdouble, (Ptr{Cbc_Model},), model)
end

function Cbc_status(model)
    ccall((:Cbc_status, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_secondaryStatus(model)
    ccall((:Cbc_secondaryStatus, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_sumPrimalInfeasibilities(model)
    ccall((:Cbc_sumPrimalInfeasibilities, libcbcsolver), Cdouble, (Ptr{Cbc_Model},), model)
end

function Cbc_numberPrimalInfeasibilities(model)
    ccall((:Cbc_numberPrimalInfeasibilities, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_checkSolution(model)
    ccall((:Cbc_checkSolution, libcbcsolver), Cvoid, (Ptr{Cbc_Model},), model)
end

function Cbc_getIterationCount(model)
    ccall((:Cbc_getIterationCount, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_isNodeLimitReached(model)
    ccall((:Cbc_isNodeLimitReached, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_isSecondsLimitReached(model)
    ccall((:Cbc_isSecondsLimitReached, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_isSolutionLimitReached(model)
    ccall((:Cbc_isSolutionLimitReached, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_isInitialSolveAbandoned(model)
    ccall((:Cbc_isInitialSolveAbandoned, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_isInitialSolveProvenOptimal(model)
    ccall((:Cbc_isInitialSolveProvenOptimal, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_isInitialSolveProvenPrimalInfeasible(model)
    ccall((:Cbc_isInitialSolveProvenPrimalInfeasible, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_getRowActivity(model)
    ccall((:Cbc_getRowActivity, libcbcsolver), Ptr{Cdouble}, (Ptr{Cbc_Model},), model)
end

function Cbc_getNodeCount(model)
    ccall((:Cbc_getNodeCount, libcbcsolver), Cint, (Ptr{Cbc_Model},), model)
end

function Cbc_printSolution(model)
    ccall((:Cbc_printSolution, libcbcsolver), Cvoid, (Ptr{Cbc_Model},), model)
end

function Osi_getNumCols(osi)
    ccall((:Osi_getNumCols, libcbcsolver), Cint, (Ptr{Cvoid},), osi)
end

function Osi_getColName(osi, i, name, maxLen)
    ccall((:Osi_getColName, libcbcsolver), Cvoid, (Ptr{Cvoid}, Cint, Cstring, Cint), osi, i, name, maxLen)
end

function Osi_getColLower(osi)
    ccall((:Osi_getColLower, libcbcsolver), Ptr{Cdouble}, (Ptr{Cvoid},), osi)
end

function Osi_getColUpper(osi)
    ccall((:Osi_getColUpper, libcbcsolver), Ptr{Cdouble}, (Ptr{Cvoid},), osi)
end

function Osi_isInteger(osi, col)
    ccall((:Osi_isInteger, libcbcsolver), Cint, (Ptr{Cvoid}, Cint), osi, col)
end

function Osi_getNumRows(osi)
    ccall((:Osi_getNumRows, libcbcsolver), Cint, (Ptr{Cvoid},), osi)
end

function Osi_getRowNz(osi, row)
    ccall((:Osi_getRowNz, libcbcsolver), Cint, (Ptr{Cvoid}, Cint), osi, row)
end

function Osi_getRowIndices(osi, row)
    ccall((:Osi_getRowIndices, libcbcsolver), Ptr{Cint}, (Ptr{Cvoid}, Cint), osi, row)
end

function Osi_getRowCoeffs(osi, row)
    ccall((:Osi_getRowCoeffs, libcbcsolver), Ptr{Cdouble}, (Ptr{Cvoid}, Cint), osi, row)
end

function Osi_getRowRHS(osi, row)
    ccall((:Osi_getRowRHS, libcbcsolver), Cdouble, (Ptr{Cvoid}, Cint), osi, row)
end

function Osi_getRowSense(osi, row)
    ccall((:Osi_getRowSense, libcbcsolver), UInt8, (Ptr{Cvoid}, Cint), osi, row)
end

function Osi_getColSolution(osi)
    ccall((:Osi_getColSolution, libcbcsolver), Ptr{Cdouble}, (Ptr{Cvoid},), osi)
end

function OsiCuts_addRowCut(osiCuts, nz, idx, coef, sense, rhs)
    ccall((:OsiCuts_addRowCut, libcbcsolver), Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}, UInt8, Cdouble), osiCuts, nz, idx, coef, sense, rhs)
end
