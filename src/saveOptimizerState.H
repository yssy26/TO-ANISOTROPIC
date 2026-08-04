//-----------------------------------------------------------------------------
// saveOptimizerState.H
// Phase J: Serialize MMA and optimization state to a file for restart.
//
// Writes to optimizerState.tmp, then atomically renames to optimizerState.dat
// to avoid corruption from interrupted writes.
//
// The saved state includes:
//   - Magic number and version for format identification
//   - Mesh/design consistency data (nCells, nDesign, activeDesignCells)
//   - MMA history (iter, xold1, xold2, low, upp, movelimit)
//   - Transaction state (currentMoveLimit, consecutiveAccepts)
//   - Continuation/projection parameters (del, eta5, alphaMax, qu)
//   - Objective and constraint values
//-----------------------------------------------------------------------------

if (mmaUpdateEnabled && enableOptimizerRestart)
{
    const fileName stateDir = runTime.path();
    const fileName tmpPath = stateDir/"optimizerState.tmp";
    const fileName datPath = stateDir/"optimizerState.dat";

    // Capture MMA state
    MMA::StateSnapshot mmaSnap = mma.saveState();

    // Capture accepted design values for active cells
    std::vector<double> xAcceptedValues(nDesignLocal);
    forAll(activeDesignCells, activeI)
    {
        xAcceptedValues[activeI] = x[activeDesignCells[activeI]];
    }

    // Write to temporary file
    OFstream os(tmpPath);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open " << tmpPath << " for writing."
            << abort(FatalError);
    }

    os.precision(16);

    // Header
    os << "MAGIC MTO_HF_MMA_STATE" << endl;
    os << "VERSION 1" << endl;

    // Mesh / design consistency
    os << "NCELLS " << n << endl;
    os << "NDESIGN_LOCAL " << nDesignLocal << endl;
    os << "NDESIGN_GLOBAL " << nDesignGlobal << endl;
    os << "NCONSTRAINTS " << m << endl;
    os << "OPT " << opt << endl;

    // Active design cells (for consistency check on load)
    os << "ACTIVE_CELLS";
    for (label i = 0; i < nDesignLocal; ++i)
    {
        os << " " << activeDesignCells[i];
    }
    os << endl;

    // MMA state
    os << "MMA_ITER " << mmaSnap.iter << endl;
    os << "MMA_MOVELIMIT " << mma.movelimit << endl;

    os << "XOLD1";
    for (label i = 0; i < nDesignLocal; ++i)
    {
        os << " " << mmaSnap.xold1[i];
    }
    os << endl;

    os << "XOLD2";
    for (label i = 0; i < nDesignLocal; ++i)
    {
        os << " " << mmaSnap.xold2[i];
    }
    os << endl;

    os << "LOW";
    for (label i = 0; i < nDesignLocal; ++i)
    {
        os << " " << mmaSnap.low[i];
    }
    os << endl;

    os << "UPP";
    for (label i = 0; i < nDesignLocal; ++i)
    {
        os << " " << mmaSnap.upp[i];
    }
    os << endl;

    // Accepted design
    os << "XACCEPTED";
    for (label i = 0; i < nDesignLocal; ++i)
    {
        os << " " << xAcceptedValues[i];
    }
    os << endl;

    // Transaction state
    os << "CURRENT_MOVE_LIMIT " << currentMoveLimit << endl;
    os << "CONSECUTIVE_ACCEPTS " << consecutiveAccepts << endl;
    os << "CSV_ACCEPTED_ITERATION " << csvAcceptedIteration << endl;

    // Continuation / projection parameters
    os << "DEL " << del << endl;
    os << "ETA5 " << eta5 << endl;
    os << "ALPHAMAX " << alphaMax.value() << endl;
    os << "QU " << qu << endl;

    // Objective and constraints
    os << "OBJECTIVE " << Objective << endl;
    os << "VOLUME_CONSTRAINT " << V << endl;
    os << "PRESSURE_DROP_PA " << PressureDropPa << endl;
    os << "MEAN_T " << MeanT << endl;
    os << "POWER_DISS " << PowerDiss << endl;

    os.flush();

    // Atomic rename
    if (std::rename(tmpPath.c_str(), datPath.c_str()) != 0)
    {
        // On some filesystems (e.g., WSL cross-OS), rename may fail.
        // Fall back to copy + remove.
        std::ifstream srcFile(tmpPath.c_str(), std::ios::binary);
        std::ofstream dstFile(datPath.c_str(), std::ios::binary);
        dstFile << srcFile.rdbuf();
        dstFile.close();
        srcFile.close();
        std::remove(tmpPath.c_str());
    }

    Info<< "Optimizer state saved to " << datPath
        << " (opt=" << opt << ", mmaIter=" << mmaSnap.iter << ")"
        << endl;
}
