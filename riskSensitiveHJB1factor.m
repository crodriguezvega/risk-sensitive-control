% riskSensitiveHJB1factor approximates solution to risk-sensitive 
% HJB equation with one factor process
%
% Reference: M.H.A. Davis and S. Lleo. Jump-diffusion risk-sensitive
% asset management I: Diffusion factor model. SIAM Journal on
% Financial Mathematics, 2:22-54, 2011.

classdef (Sealed = true) riskSensitiveHJB1factor
properties (GetAccess = private, SetAccess = private)
  % Value function
  V
  h
  % Wealth dynamics
  a_tilde
  A_tilde
  % Control policy
  weights
  isControlAdmissible
  % Results from integration with respect to Levy measure
  integrals
end
properties
  t        % Time range
  x1       % Factor value range
  numberOfControls
  % Risk sensitivity
  theta    
  % Factor dynamics 
  b
  B
  lambda
  % Asset market dynamics
  a0
  A0 
  % Risky security dynamics
  a
  A
  sigma
  z
  gamma_min
  gamma_max
  gamma
  % Small/big jump border
  R
  % HJB problem 
  v
  % Levy measure with Gaussian density
  intensity
  mean
  variance
end

methods
function obj = riskSensitiveHJB1factor()
end

function [V, h] = solve(obj)  
% Check if input values are correct
obj.checkParameters();

dt = diff(obj.t(1:2));
dx1 = diff(obj.x1(1:2));
    
obj.a_tilde = obj.a - obj.a0;
obj.A_tilde = obj.A - obj.A0;

% Preallocate matrices
% Value function
obj.V = zeros(length(obj.x1) + 2, length(obj.t)); 
% Optimal control
obj.h = zeros(length(obj.x1) + 2, length(obj.t) - 1); 

% Calculate asset weights
obj.weights = linspace(-1/obj.gamma_max, 1, obj.numberOfControls + 1);
obj.weights = obj.weights(2:end);

% Find what controls are admissible
obj.isControlAdmissible = obj.findAdmissibleControls();

% Precompute integrals with respect to Levy measure
obj.integrals = obj.calculateIntegralsWithRespectToLevyMeasure();

% Apply terminal condition to value function
obj.V(:, length(obj.t)) = log(obj.v);

for m = length(obj.t):-1:2  % time
  if m < length(obj.t)
    obj.V(:, m) = obj.extrapolateValueFunctionBeyondBorders(dx1, m);
  end
  
  for i = 2:1:length(obj.x1) + 1  % factor
    [DV, D2V] = obj.partialDerivatives(dx1, m, i);    
    [sup, optimal_h] = obj.supOperatorL(obj.x1(i - 1), DV); 
    obj.h(i, m - 1, :) = optimal_h;
    obj.V(i, m - 1) = obj.V(i, m) ...
                      + dt*((obj.b + obj.B*obj.x1(i - 1))*(DV) ...
                      + 0.5*(obj.lambda^2)*(D2V) ...
                      - (obj.theta/2)*(obj.lambda^2)*(DV^2) ...
                      + obj.a0 + obj.A0*obj.x1(i - 1) + sup);    

    fprintf('Value function for state %.2f calculated at t = %.4f\n', obj.x1(i - 1), obj.t(m - 1));
  end
end
V = obj.V(2:end - 1, :);
h = obj.h(2:end - 1, :, :);
end
end % public methods

methods (Access = private)
%CHECKPARAMETERS Check that some of the input parameters are valid
% in order to solve the HJB equation.
function checkParameters(obj)
if obj.gamma_min <= -1 || obj.gamma_min >= 0
  error('riskSensitiveHJB1factor:invalidInputs', 'gamma_min must be > -1 and < 0.');
end
if obj.gamma_max <= 0
  error('riskSensitiveHJB1factor:invalidInputs', 'gamma_max must be > 0.');
end
if obj.theta == 0 || obj.theta <= -1
  error('riskSensitiveHJB1factor:invalidInputs', 'theta cannot be 0 or <= -1.');
end
end

%EXTRAPOLATEVALUEFUNCTIONBEYONDBORDERS Extrapolates one extra value 
% around the borders of the value function
function [val] = extrapolateValueFunctionBeyondBorders(obj, dx1, m)
Xq = obj.x1(1) - dx1:dx1:obj.x1(end) + dx1;
val = interp1(obj.x1, obj.V(2:end - 1, m), Xq, 'linear', 'extrap');
end
    
%PARTIALDERIVATIVES Return first and second order partial derivatives
% with respect to the state variable
function [DV, D2V] = partialDerivatives(obj, dx1, m, i)
% First order partial derivative of value function
DV = (obj.V(i + 1, m) - obj.V(i - 1, m))/(2*dx1); 

% Second order partial derivative of value function  
D2V = (obj.V(i + 1, m) - 2*obj.V(i, m) + obj.V(i - 1, m))/(dx1^2);
end

%SUPOPERATORL Calculate the maximum value and the control policy that
% produces it
function [val, optimal_h] = supOperatorL(obj, x, DV)
f = zeros(length(obj.weights), 1);           
for i = 1:1:length(obj.weights)   % Weight for risky asset
  if obj.isControlAdmissible(i) == true
    h = obj.weights(i);

    f(i) = -0.5*(obj.theta + 1)*(h^2)*(obj.sigma^2) ...
           - obj.theta*(h*obj.sigma*obj.lambda*DV) ...
           + h*(obj.a_tilde + obj.A_tilde*x) ...
           - (1/obj.theta)*obj.integrals(i);               
  else
    f(i) = nan;
  end
end

% Find maximum value
[val, ind] = max(f(:));
[row, column] = ind2sub(size(f), ind); 
% Control policy that produces the maximum
optimal_h = obj.weights(row);
end

%FINDADMISSIBLECONTROLS Find what control policies are admissible
function [isControlAdmissible] = findAdmissibleControls(obj)
isControlAdmissible = zeros(length(obj.weights), 1);           
for i = 1:1:length(obj.weights)   % Weight for risky asset
  if obj.checkControlIsAdmissible(obj.weights(i)) == true                
    isControlAdmissible(i) = true;
  else
    isControlAdmissible(i) = false;
  end         
end
end

%CHECKCONTROLISADMISSIBLE Check if input control policy is admissible
function [isAdmissible] = checkControlIsAdmissible(obj, h)
isAdmissible = true;
for i = 1:1:length(obj.gamma)
  % Condition for admissibility
  if obj.gamma(i)*h <= -1 || h > 1
    isAdmissible = false;
    break;
  end                
end
end

%CALCULATEINTEGRALWITHRESPECTTOLEVYMEASURE Calculate integral with
% respect to Levy measure for every admissible control policy
function [integrals] = calculateIntegralsWithRespectToLevyMeasure(obj)
fprintf('Calculating integrals with respect to Levy measure');
  
integrals = zeros(length(obj.weights), 1);       
for i = 1:1:length(obj.weights)   % Weight for risky asset  
  if obj.isControlAdmissible(i) == true   
    integrals(i) = obj.integralWithRespectToLevyMeasure(obj.weights(i));
  else
    integrals(i) = nan;
  end
  fprintf('.');
end

fprintf('\n');
end

%INTEGRALWITHRESPECTTOLEVYMEASURE Approximate integral with respect to
% Levy measure
function val = integralWithRespectToLevyMeasure(obj, h)
val = 0;
dz = diff(obj.z(1:2));

for i = 2:1:length(obj.z) - 1
  % Evaluate the integrand and accumulate result  
  val = val + obj.evaluateIntegrand(i, h);     
end
val = dz*(obj.evaluateIntegrand(1, h)/2 + val + obj.evaluateIntegrand(length(obj.z), h)/2);
end

%EVALUATEINTEGRAND Evaluate integrand of integral with respect to
% Levy measure at one point of the domain
function val = evaluateIntegrand(obj, i, h)
aux1 = (power(1 + obj.gamma(i)*h, -1*obj.theta) - 1);
aux2 = 0;
if abs(obj.z(i)) <= obj.R
  aux2 = obj.theta*obj.gamma(i)*h;
end
val = (aux1 + aux2)*obj.gaussianDensity1D(obj.z(i));
end

%GAUSSIANDENSITY1D Calculate the value of the 1D Gaussian density
% function at any input point
function val = gaussianDensity1D(obj, z)
if z == 0
  val = 0; % Levy measure has no mass at origin
else
  val = obj.intensity*normpdf(z, obj.mean, obj.variance);
end
end

end % private methods
end % classdef