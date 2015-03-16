function tests = test_nlcf
    tests = functiontests({
        @testInputValidation
        @testContinuousTimeAlgorithm
        @testDiscreteTimeAlgorithm
    });
end

function testInputValidation(testCase)
% Test input validation
    import cftool.nlcf;
    testCase.verifyError(@nlcf, 'MATLAB:narginchk:notEnoughInputs');
    testCase.verifyError(@() nlcf([]), 'MATLAB:invalidType');
end

function testDiscreteTimeAlgorithm(testCase)
% Test the discrete-time normalized left coprime factorization.
    runTest(testCase, @drss);
end

function testContinuousTimeAlgorithm(testCase)
% Test the continuous time normalized left coprime factorization.
    runTest(testCase, @rss);
end

function runTest(testCase, sysGen)
    dropIn = true;
    function run(n, nu, ny)
        import cftool.*;

        % Run over this many random systems.
        Nruns = 8;

        % Continuous time
        for i = 1:Nruns
            G0 = sysGen(n, ny, nu);
            try
                [N, M] = nlcf(G0);
            catch err
                if dropIn
                    keyboard;
                else
                    rethrow(err)
                end
            end
            
            if ~isIoEqual(M\N, G0)
                keyboard;
                testCase.assertFail;
            end
            
            R = minreal(balreal(M*M' + N*N'), [], false);
            testCase.assertLessThan(abs(dcgain(R) - eye(ny)), 1e-7);
        end
    end

    for n = [1, 2, 6]
        for nu = [1, 2, 6]
            for ny = [1, 2, 6]
                run(n, nu, ny);
            end
        end
    end
end
