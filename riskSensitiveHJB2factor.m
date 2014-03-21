% riskSensitiveHJB2factor approximates solution to risk-sensitive
% HJB equation with two factor process
%
% Reference: M.H.A. Davis and S. Lleo. Jump-diffusion risk-sensitive
% asset management I: Diffusion factor model. SIAM Journal on
% Financial Mathematics, 2:22-54, 2011.

classdef (Sealed = true) riskSensitiveHJB2factor
properties (Constant)
  one = [1
         1];
end
properties (GetAccess = private, SetAccess = private)
  % Value function
  V
  % Optimal control
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
  x1       % Factor 1 value range
  x2       % Factor 2 value range
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
  z1
  z2
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
  covariance
end

methods
function obj = riskSensitiveHJB2factor()
end

function [V, h] = solve(obj)
% Check if input values are correct
obj.checkParameters();

dt = diff(obj.t(1:2));
dx1 = diff(obj.x1(1:2));
dx2 = diff(obj.x2(1:2));
    
obj.a_tilde = obj.a - obj.a0.*riskSensitiveHJB2factor.one;
obj.A_tilde = obj.A - riskSensitiveHJB2factor.one*obj.A0';

% Preallocate matrices
% Value function
obj.V = zeros(length(obj.x1) + 2, length(obj.x2) + 2, length(obj.t));
% Optimal control
obj.h = zeros(length(obj.x1) + 2, length(obj.x2) + 2, length(obj.t) - 1, 2);

% Calculate asset weights
obj.weights = linspace(-1/obj.gamma_max, 1, obj.numberOfControls + 1);
obj.weights = obj.weights(2:end);

% Find what controls are admissible
obj.isControlAdmissible = obj.findAdmissibleControls();

% Precompute integrals with respect to Levy measure
obj.integrals = obj.calculateIntegralsWithRespectToLevyMeasure();

% Apply terminal condition to value function
obj.V(:, :, length(obj.t)) = log(obj.v);

for m = length(obj.t):-1:2          % time
  if m < length(obj.t)
    obj.V(:, :, m) = obj.extrapolateValueFunctionBeyondBorders(dx1, dx2, m);
  end
  
  for i = 2:1:length(obj.x1) + 1    % factor 1
    for j = 2:1:length(obj.x2) + 1  % factor 2
      x = [obj.x1(i - 1)
           obj.x2(j - 1)];

      [DV, D2V] = obj.partialDerivatives(dx1, dx2, m, i, j);
      [sup, optimal_h] = obj.supOperatorL(x, DV); 
      obj.h(i, j, m - 1, :) = optimal_h;
      obj.V(i, j, m - 1) = obj.V(i, j, m) ...
        + dt*((obj.b + obj.B*x)'*(DV) ...
        + 0.5*trace(obj.lambda*(obj.lambda'*D2V)) ...
        - (obj.theta/2)*(DV'*(obj.lambda*(obj.lambda'*DV))) ...
        + (obj.a0 + obj.A0'*x + sup));    

      fprintf('Value function for state (%.2f, %.2f) calculated at t = %.4f\n', obj.x1(i - 1), obj.x2(j - 1), obj.t(m - 1));
    end
  end
end
V = obj.V(2:end - 1, 2:end - 1, :);
h = obj.h(2:end - 1, 2:end - 1, :, :);
end
end % public methods

methods (Access = private)
%CHECKPARAMETERS Check that some of the input parameters are valid 
% in order to solve the HJB equation.
function checkParameters(obj)
% Check that lambda*lambda' is positive definite
[R, p] = chol(obj.lambda*obj.lambda');
if p > 0
  error('riskSensitiveHJB2factor:invalidInputs', 'lambda*lambda'' is not positive definite.');
end
% Check that sigma*sigma' is positive definite
[R, p] = chol(obj.sigma*obj.sigma');
if p > 0
  error('riskSensitiveHJB2factor:invalidInputs', 'sigma*sigma'' is not positive definite.');
end
if obj.gamma_min <= -1 || obj.gamma_min >= 0
  error('riskSensitiveHJB2factor:invalidInputs', 'gamma_min must be > -1 and < 0.');
end
if obj.gamma_max <= 0
  error('riskSensitiveHJB2factor:invalidInputs', 'gamma_max must be > 0.');
end
if obj.theta == 0 || obj.theta <= -1
  error('riskSensitiveHJB2factor:invalidInputs', 'theta cannot be 0 or <= -1.');
end
end
    
%EXTRAPOLATEVALUEFUNCTIONBEYONDBORDERS Extrapolates one extra value 
% around the borders of the value function
function [val] = extrapolateValueFunctionBeyondBorders(obj, dx1, dx2, m)
[X, Y] = meshgrid(obj.x2, obj.x1);
[Xq, Yq] = meshgrid(obj.x2(1) - dx2:dx2:obj.x2(end) + dx2, obj.x1(1) - dx1:dx1:obj.x1(end) + dx1);  
val = interp2(X, Y, obj.V(2:end - 1, 2:end - 1, m), Xq, Yq, 'spline');
end

%PARTIALDERIVATIVES Return first and second order partial derivatives
% with respect to the state variable
function [DV, D2V] = partialDerivatives(obj, dx1, dx2, m, i, j)
% First order partial derivative of value function
DV = [(obj.V(i + 1, j, m) - obj.V(i - 1, j, m))/(2*dx1) 
      (obj.V(i, j + 1, m) - obj.V(i, j - 1, m))/(2*dx2)]; 

% Second order partial derivative of value function (here x=x1, y=x2)  
fxx = (obj.V(i + 1, j, m) - 2*obj.V(i, j, m) ...
  + obj.V(i - 1, j, m))/(dx1^2);
fyy = (obj.V(i, j + 1, m) - 2*obj.V(i, j, m) ...
  + obj.V(i, j - 1, m))/(dx2^2);
fxy = (obj.V(i + 1, j + 1, m) - obj.V(i + 1, j - 1, m) ...
  - obj.V(i - 1, j + 1, m) + obj.V(i - 1, j - 1, m))/(4*dx1*dx2);
fyx = fxy;            
D2V = [fxx fxy
       fyx fyy];
end

%FINDADMISSIBLECONTROLS Find what control policies are admissible
function [isControlAdmissible] = findAdmissibleControls(obj)
isControlAdmissible = zeros(length(obj.weights), length(obj.weights));           
for i = 1:1:length(obj.weights)   % Weight for asset 1
  for j = 1:1:length(obj.weights) % Weight for asset 2
    h = [obj.weights(i)
         obj.weights(j)];

    if obj.checkControlIsAdmissible(h) == true                
      isControlAdmissible(i, j) = true;
    else
      isControlAdmissible(i, j) = false;
    end         
  end
end
end

%CHECKCONTROLISADMISSIBLE Check if input control policy is admissible
function [isAdmissible] = checkControlIsAdmissible(obj, h)
isAdmissible = true;
for i = 1:1:length(obj.gamma(:, 1, 1))
  if isAdmissible == false
    break;
  end

  for j = 1:1:length(obj.gamma(1, :, 1))
    % Condition for admissibility
    if [obj.gamma(i, j, 1) obj.gamma(i, j, 2)]*h <= -1 || sum(h) > 1
      isAdmissible = false;
      break;
    end                
  end
end
end

%SUPOPERATORL Calculate the maximum value and the control policy 
% that produces it
function [val, optimal_h] = supOperatorL(obj, x, DV)
f = zeros(length(obj.weights), length(obj.weights));           
for i = 1:1:length(obj.weights)   % Weight for asset 1
  for j = 1:1:length(obj.weights) % Weight for asset 2    
    if obj.isControlAdmissible(i, j) == true
      h = [obj.weights(i)
           obj.weights(j)];

      f(i, j) = -0.5*(obj.theta + 1)*(h'*obj.sigma)*(obj.sigma'*h) ...
                - obj.theta*(h'*(obj.sigma*(obj.lambda'*DV))) ...
                + h'*(obj.a_tilde + obj.A_tilde*x) ...
                - (1/obj.theta)*obj.integrals(i, j);               
    else
      f(i, j) = nan;
    end
  end
end

% Find maximum value
[val, ind] = max(f(:));
[row, column] = ind2sub(size(f), ind); 
% Control policy that produces the maximum
optimal_h = [obj.weights(row)
             obj.weights(column)];
end

%CALCULATEINTEGRALWITHRESPECTTOLEVYMEASURE Calculate integral with
% respect to Levy measure for every admissible control policy
function [integrals] = calculateIntegralsWithRespectToLevyMeasure(obj)
fprintf('Calculating integrals with respect to Levy measure');

integrals = zeros(length(obj.weights), length(obj.weights));       
for i = 1:1:length(obj.weights)   % Weight for asset 1  
  for j = 1:1:length(obj.weights) % Weight for asset 2
    if obj.isControlAdmissible(i, j) == true
      h = [obj.weights(i)
           obj.weights(j)];

      integrals(i, j) = obj.integralWithRespectToLevyMeasure(h);
    else
      integrals(i, j) = nan;
    end
  end
  fprintf('.');
end

fprintf('\n');
end

%INTEGRALWITHRESPECTTOLEVYMEASURE Approximate integral with respect
% to Levy measure
function val = integralWithRespectToLevyMeasure(obj, h)
val = 0;
dz1 = diff(obj.z1(1:2));
dz2 = diff(obj.z2(1:2));

for i = 1:1:length(obj.z1) - 1
  for j = 1:1:length(obj.z2) - 1
    % Evaluate the integrand at every corner of the rectangle  
    corner1 = obj.evaluateIntegrand(i, j, h);
    corner2 = obj.evaluateIntegrand(i + 1, j, h);
    corner3 = obj.evaluateIntegrand(i, j + 1, h);
    corner4 = obj.evaluateIntegrand(i + 1, j + 1, h);
    val = val + 0.25*(corner1 + corner2 + corner3 + corner4)*dz1*dz2;     
  end
end
end

%EVALUATEINTEGRAND Evaluate integrand of integral with respect to
% Levy measure at one point of the domain
function val = evaluateIntegrand(obj, i, j, h)
aux1 = (power(1 + [obj.gamma(i, j, 1) obj.gamma(i, j, 2)]*h, -1*obj.theta) - 1);
aux2 = 0;
if (abs(obj.z1(i)) <= obj.R(1) && abs(obj.z2(j)) <= obj.R(2))
  aux2 = obj.theta*[obj.gamma(i, j, 1) obj.gamma(i, j, 2)]*h;
end
z = [obj.z1(i)
     obj.z2(j)];
val = (aux1 + aux2)*obj.gaussianDensity2D(z);
end

%GAUSSIANDENSITY2D Calculate the value of the 2D Gaussian density
% function at any input point
function val = gaussianDensity2D(obj, z)
if isequal(z, [0 0]')
  val = 0; % Levy measure has no mass at origin
else
  val = obj.intensity*(1/sqrt(power(2*pi, length(obj.mean))*det(obj.covariance)))*(exp(-0.5*(z - obj.mean)'*(obj.covariance\(z - obj.mean))));
end
end

end % private methods
end % classdef