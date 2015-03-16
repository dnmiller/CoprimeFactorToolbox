function tests = test_nrcf
    tests = functiontests({
        @testInputValidation
        @testContinuousTimeAlgorithm
        @testDiscreteTimeAlgorithm
    });
end

function testInputValidation(testCase)
% Test input validation
    import cftool.nrcf;
    testCase.verifyError(@nrcf, 'MATLAB:narginchk:notEnoughInputs');
    testCase.verifyError(@() nrcf([]), 'MATLAB:invalidType');
end

function testDiscreteTimeAlgorithm(testCase)
    dropIn = false;
    function run_test(n, nu, ny)
        import cftool.*;

        % Run over this many random systems.
        Nruns = 50;

        % Continuous time
        for i = 1:Nruns
            G0 = drss(n, ny, nu);
            try
                [N, M] = nrcf(G0);
            catch err
                if dropIn
                    keyboard;
                    [N, M] = nrcf(G0);
                else
                    rethrow(err)
                end
            end
            
            if ~isIoEqual(N/M, G0)
                keyboard;
                testCase.assertFail;
            end
            
            R = minreal(balreal(M'*M + N'*N), [], false);
            testCase.assertLessThan(abs(dcgain(R) - eye(nu)), 1e-7);
        end
    end

    for n = [1, 2, 6]
        for nu = [1, 2, 6]
            for ny = [1, 2, 6]
                run_test(n, nu, ny);
            end
        end
    end
end

function testContinuousTimeAlgorithm(testCase)
% Test the continuous time normalized right coprime factorization.
return
    dropIn = false;
    function run_test(n, nu, ny)
        import cftool.*;

        % Run over this many random systems.
        Nruns = 8;

        % Continuous time
        for i = 1:Nruns
            G0 = rss(n, ny, nu);
            try
                [N, M] = nrcf(G0);
            catch err
                if dropIn
                    keyboard;
                    [N, M] = nrcf(G0);
                else
                    rethrow(err)
                end
            end
            
            if ~isIoEqual(N/M, G0)
                keyboard;
                testCase.assertFail;
            end
            
            R = minreal(balreal(M'*M + N'*N), [], false);
            testCase.assertLessThan(abs(dcgain(R) - eye(nu)), sqrt(eps));
        end
    end

    for n = [1, 2, 6]
        for nu = [1, 2, 6]
            for ny = [1, 2, 6]
                run_test(n, nu, ny);
            end
        end
    end
end