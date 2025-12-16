classdef bfUnitTest < matlab.unittest.TestCase
    % Class to test the Bayes Factor toolbox using the Matlab Unit Test
    % framework.
    % BK - Dec 2019
    
    methods (TestClassSetup)
        function importPackages(testCase) %#ok<MANU>
            import matlab.unittest.*
            import matlab.unittest.constraints.*            
        end
    end
    
    methods (Test)
        function testFromF(testCase)
            actual = bf.bfFromF(10,1,3,5);
            expected = 17.4812;
            testCase.verifyThat(actual,matlab.unittest.constraints.IsEqualTo(expected,'Within',matlab.unittest.constraints.AbsoluteTolerance(0.0001)))
        end
        
        function testFromT(testCase)
            actual = bf.bfFromT(2,5);
            expected = 2.9574;
            testCase.verifyThat(actual,matlab.unittest.constraints.IsEqualTo(expected,'Within',matlab.unittest.constraints.AbsoluteTolerance(0.0001)))
        end
        
        function testTTestPValueTwoTailed(testCase)
            % Test the bug fix: two-tailed test with T=3.04, N=30
            % should produce p-value < 1 (specifically around 0.005)
            [~, pValue] = bf.ttest('T', 3.04, 'N', 30, 'tail', 'both');
            testCase.verifyLessThan(pValue, 1, 'P-value should be less than 1');
            testCase.verifyGreaterThan(pValue, 0, 'P-value should be greater than 0');
            % For T=3.04, df=29, two-tailed p should be approximately 0.0050
            testCase.verifyThat(pValue, matlab.unittest.constraints.IsEqualTo(0.0050, 'Within', matlab.unittest.constraints.AbsoluteTolerance(0.001)));
        end
        
        function testTTestPValueRightTailed(testCase)
            % Test right-tailed p-value calculation
            [~, pValue] = bf.ttest('T', 2.5, 'N', 20, 'tail', 'right');
            testCase.verifyLessThan(pValue, 1, 'P-value should be less than 1');
            testCase.verifyGreaterThan(pValue, 0, 'P-value should be greater than 0');
            % For positive T, right-tailed should give small p-value
            testCase.verifyLessThan(pValue, 0.05, 'P-value should be small for positive T in right-tailed test');
        end
        
        function testTTestPValueLeftTailed(testCase)
            % Test left-tailed p-value calculation
            [~, pValue] = bf.ttest('T', -2.5, 'N', 20, 'tail', 'left');
            testCase.verifyLessThan(pValue, 1, 'P-value should be less than 1');
            testCase.verifyGreaterThan(pValue, 0, 'P-value should be greater than 0');
            % For negative T, left-tailed should give small p-value
            testCase.verifyLessThan(pValue, 0.05, 'P-value should be small for negative T in left-tailed test');
        end
        
        function testTTestPValueNegativeTTwoTailed(testCase)
            % Test two-tailed test with negative T value
            [~, pValue] = bf.ttest('T', -3.04, 'N', 30, 'tail', 'both');
            testCase.verifyLessThan(pValue, 1, 'P-value should be less than 1');
            testCase.verifyGreaterThan(pValue, 0, 'P-value should be greater than 0');
            % Should be approximately same as positive T for two-tailed
            testCase.verifyThat(pValue, matlab.unittest.constraints.IsEqualTo(0.0050, 'Within', matlab.unittest.constraints.AbsoluteTolerance(0.001)));
        end

    end
end