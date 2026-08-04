//-----------------------------------------------------------------------------
// loadOptimizerState.H
// Phase J: Deserialize MMA and optimization state from optimizerState.dat.
//
// Performs strict consistency checks:
//   - Magic number and version
//   - nCells, nDesignLocal, nDesignGlobal, m
//   - Active design cells match current mesh
//
// On any inconsistency, the state file is rejected and the optimization
// starts from scratch.
//
// This file is included from opt_initialization.H after the MMA object is
// created.  It is only active when enableOptimizerRestart is true.
//-----------------------------------------------------------------------------

if (mmaUpdateEnabled && enableOptimizerRestart)
{
    const fileName datPath = runTime.path()/"optimizerState.dat";

    if (isFile(datPath))
    {
        Info<< "Loading optimizer state from " << datPath << endl;

        IFstream is(datPath);
        if (!is.good())
        {
            FatalErrorInFunction
                << "Cannot open " << datPath << " for reading."
                << abort(FatalError);
        }

        // --- Magic and version ---
        word magicWord, magicValue;
        is >> magicWord >> magicValue;

        if (magicWord != "MAGIC" || magicValue != "MTO_HF_MMA_STATE")
        {
            FatalErrorInFunction
                << "Invalid magic number in " << datPath
                << ": expected 'MTO_HF_MMA_STATE', got '" << magicValue
                << "'. State file rejected." << abort(FatalError);
        }

        word versionWord;
        label version;
        is >> versionWord >> version;

        if (versionWord != "VERSION" || version != 1)
        {
            FatalErrorInFunction
                << "Unsupported state file version " << version
                << ". Expected version 1." << abort(FatalError);
        }

        // --- Mesh / design consistency ---
        word key;
        label fileNCells, fileNDesignLocal, fileNDesignGlobal, fileM, fileOpt;

        is >> key >> fileNCells;
        if (key != "NCELLS" || fileNCells != n)
        {
            FatalErrorInFunction
                << "nCells mismatch: file=" << fileNCells
                << ", mesh=" << n << abort(FatalError);
        }

        is >> key >> fileNDesignLocal;
        if (key != "NDESIGN_LOCAL" || fileNDesignLocal != nDesignLocal)
        {
            FatalErrorInFunction
                << "nDesignLocal mismatch: file=" << fileNDesignLocal
                << ", current=" << nDesignLocal << abort(FatalError);
        }

        is >> key >> fileNDesignGlobal;
        if (key != "NDESIGN_GLOBAL" || fileNDesignGlobal != nDesignGlobal)
        {
            FatalErrorInFunction
                << "nDesignGlobal mismatch: file=" << fileNDesignGlobal
                << ", current=" << nDesignGlobal << abort(FatalError);
        }

        is >> key >> fileM;
        if (key != "NCONSTRAINTS" || fileM != m)
        {
            FatalErrorInFunction
                << "nConstraints mismatch: file=" << fileM
                << ", current=" << m << abort(FatalError);
        }

        is >> key >> fileOpt;
        if (key != "OPT")
        {
            FatalErrorInFunction
                << "Expected OPT key, got " << key << abort(FatalError);
        }

        // --- Active design cells ---
        is >> key;
        if (key != "ACTIVE_CELLS")
        {
            FatalErrorInFunction
                << "Expected ACTIVE_CELLS key, got " << key
                << abort(FatalError);
        }

        bool activeCellsMatch = true;
        for (label i = 0; i < nDesignLocal; ++i)
        {
            label fileCell;
            is >> fileCell;
            if (fileCell != activeDesignCells[i])
            {
                activeCellsMatch = false;
            }
        }

        if (!activeCellsMatch)
        {
            FatalErrorInFunction
                << "Active design cells do not match the current mesh."
                << " The state file was generated for a different mesh"
                << " or design mask configuration." << abort(FatalError);
        }

        // --- MMA state ---
        MMA::StateSnapshot mmaSnap;
        mmaSnap.xold1.resize(nDesignLocal);
        mmaSnap.xold2.resize(nDesignLocal);
        mmaSnap.low.resize(nDesignLocal);
        mmaSnap.upp.resize(nDesignLocal);

        is >> key >> mmaSnap.iter;
        if (key != "MMA_ITER")
        {
            FatalErrorInFunction
                << "Expected MMA_ITER, got " << key << abort(FatalError);
        }

        scalar fileMmaMove;
        is >> key >> fileMmaMove;
        if (key != "MMA_MOVELIMIT")
        {
            FatalErrorInFunction
                << "Expected MMA_MOVELIMIT, got " << key << abort(FatalError);
        }

        // Read xold1, xold2, low, upp
        auto readVector = [&](word expectedKey, std::vector<double>& vec)
        {
            is >> key;
            if (key != expectedKey)
            {
                FatalErrorInFunction
                    << "Expected " << expectedKey << ", got " << key
                    << abort(FatalError);
            }
            for (label i = 0; i < nDesignLocal; ++i)
            {
                is >> vec[i];
            }
        };

        readVector("XOLD1", mmaSnap.xold1);
        readVector("XOLD2", mmaSnap.xold2);
        readVector("LOW", mmaSnap.low);
        readVector("UPP", mmaSnap.upp);

        // Read accepted design
        std::vector<double> fileXAccepted(nDesignLocal);
        readVector("XACCEPTED", fileXAccepted);

        // --- Transaction state ---
        scalar fileMoveLimit;
        label fileConsecutiveAccepts;

        is >> key >> fileMoveLimit;
        if (key != "CURRENT_MOVE_LIMIT")
        {
            FatalErrorInFunction
                << "Expected CURRENT_MOVE_LIMIT, got " << key
                << abort(FatalError);
        }

        is >> key >> fileConsecutiveAccepts;
        if (key != "CONSECUTIVE_ACCEPTS")
        {
            FatalErrorInFunction
                << "Expected CONSECUTIVE_ACCEPTS, got " << key
                << abort(FatalError);
        }

        label fileCsvAcceptedIter;
        is >> key >> fileCsvAcceptedIter;
        if (key != "CSV_ACCEPTED_ITERATION")
        {
            FatalErrorInFunction
                << "Expected CSV_ACCEPTED_ITERATION, got " << key
                << abort(FatalError);
        }

        // --- Continuation / projection parameters ---
        auto checkKey = [&](const word& expected, const word& actual)
        {
            if (actual != expected)
            {
                FatalErrorInFunction
                    << "Expected " << expected << ", got " << actual
                    << " in optimizer state file." << abort(FatalError);
            }
        };

        scalar fileDel, fileEta5, fileAlphaMax, fileQu;

        is >> key >> fileDel;
        checkKey("DEL", key);

        is >> key >> fileEta5;
        checkKey("ETA5", key);

        is >> key >> fileAlphaMax;
        checkKey("ALPHAMAX", key);

        is >> key >> fileQu;
        checkKey("QU", key);

        // --- Objective and constraints ---
        scalar fileObjective, fileV, filePressureDrop, fileMeanT, filePowerDiss;

        is >> key >> fileObjective;
        checkKey("OBJECTIVE", key);

        is >> key >> fileV;
        checkKey("VOLUME_CONSTRAINT", key);

        is >> key >> filePressureDrop;
        checkKey("PRESSURE_DROP_PA", key);

        is >> key >> fileMeanT;
        checkKey("MEAN_T", key);

        is >> key >> filePowerDiss;
        checkKey("POWER_DISS", key);

        // --- Apply loaded state ---
        mma.restoreState(mmaSnap);
        mma.movelimit = fileMmaMove;

        currentMoveLimit = fileMoveLimit;
        consecutiveAccepts = fileConsecutiveAccepts;
        csvAcceptedIteration = fileCsvAcceptedIter;

        del = fileDel;
        eta5 = fileEta5;
        alphaMax.value() = fileAlphaMax;
        qu = fileQu;

        Objective = fileObjective;
        V = fileV;
        PressureDropPa = filePressureDrop;
        MeanT = fileMeanT;
        PowerDiss = filePowerDiss;

        // Verify accepted design matches the loaded x field
        scalar maxDesignDiff = 0.0;
        forAll(activeDesignCells, activeI)
        {
            const scalar diff = mag
            (
                x[activeDesignCells[activeI]] - fileXAccepted[activeI]
            );
            maxDesignDiff = Foam::max(maxDesignDiff, diff);
        }
        reduce(maxDesignDiff, maxOp<scalar>());

        if (maxDesignDiff > 1e-6)
        {
            WarningInFunction
                << "Design field differs from saved optimizer state by "
                << maxDesignDiff << " (max). This may indicate the restart"
                << " time step does not match the saved state." << endl;
        }

        Info<< "Optimizer state loaded successfully:" << endl;
        Info<< "  opt=" << fileOpt << endl;
        Info<< "  mmaIter=" << mmaSnap.iter << endl;
        Info<< "  currentMoveLimit=" << currentMoveLimit << endl;
        Info<< "  consecutiveAccepts=" << consecutiveAccepts << endl;
        Info<< "  csvAcceptedIteration=" << csvAcceptedIteration << endl;
        Info<< "  objective=" << Objective << endl;
        Info<< "  maxDesignDiff=" << maxDesignDiff << endl;
    }
    else
    {
        Info<< "No optimizer state file found at " << datPath
            << ". Starting optimization from scratch." << endl;
    }
}
