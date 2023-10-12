classdef sysall < matlab.System
    
    % System object implementation of a 2-state EKF for attitude estimation
    % For more info on the algorithm, its implementation and usage, go to
    % the root folder and see the files README.txt, Presentation.pdf, and
    % ekf2sim_example.slx (read the help of the ekf2sim block therein).
    % Copyright 2014 The MathWorks, Inc.
    
    properties
        W =     1.e-5*[0.1134 0; 0 0.3123];
        V =     [0.3368 0; 0 0.092];
        P0 =    zeros(2);
        x0 =    [0.1667 0.3825]';
        Plim =  1e-5;
        T =     0.05;
        V0 =    [-21.4287 -0.7476 8.0028]';
    end
    
    properties (DiscreteState)
        P
        x
        vold
    end
    
    methods
        function obj = sysall(varargin)
            % Support name-value pair arguments
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access=protected)
        
        function resetImpl(obj)
            % Specify initial values for DiscreteState properties
            obj.P=obj.P0;
            obj.x=obj.x0(:);
            obj.vold=obj.V0(:);
        end
        
        function rpy = stepImpl(obj,u,ab,ve)
            % Implement System algorithm:
            
            %% Calculate rpy as a function of inputs and state.
            
            % calculate heading
            psi=atan2(ve(2),ve(1));
            
            % attach psi to x and send all to output
            rpy = [obj.x;psi];
            
            
            %% Calculate next state
            
            % calculate gps acceleration
            agps = (ve-obj.vold)/obj.T;
            
            % rotate gps acceleration
            rx=-cos(psi)*agps(1)-sin(psi)*agps(2);
            ry=sin(psi)*agps(1)-cos(psi)*agps(2);
            rz=9.80665-agps(3);
            
            % solve for theta and phi given imu and gps accelerations
            sth=(rx*ab(1)+rz*(abs(rx^2+rz^2-ab(1)^2)^0.5))/(realmin+rx^2+rz^2);
            theta=atan2(sth*rx-ab(1),sth*rz);
            rth=rx*sin(theta)+rz*cos(theta);
            sph=(ry*ab(2)+rth*(abs(ry^2+rth^2-ab(2)^2)^0.5))/(realmin+ry^2+rth^2);
            phi=atan2(-sph*ry+ab(2),sph*rth);
            
            % this works as the "measured" output in the ekf
            y = [phi;theta];
            
            % states (note, phi and theta are re-initialized with the previous estimates)
            phi=obj.x(1);
            theta=obj.x(2);
            
            % inputs
            p=u(1);
            q=u(2);
            r=u(3);
            
            % omega to attitude dot
            Taw=[1 sin(phi)*tan(theta) cos(phi)*tan(theta); 0 cos(phi) -sin(phi)];
            
            % x propagation
            xdot=Taw*[p;q;r];
            xm=obj.x+obj.T*xdot;
            
            % A matrix
            a = [ cos(phi)*tan(theta)*q-sin(phi)*tan(theta)*r      sin(phi)*(1+tan(theta)^2)*q+cos(phi)*(1+tan(theta)^2)*r
                -sin(phi)*q-cos(phi)*r                           0                                                      ];
            
            % euler forward step
            A = eye(2)+obj.T*a;
            
            % find K
            Pm=obj.W+A*obj.P*A';
            VPm=obj.V+Pm;
            VPm_det=VPm(1,1)*VPm(2,2)-VPm(1,2)*VPm(2,1);
            iVPm=[VPm(2,2) -VPm(1,2); -VPm(2,1) VPm(1,1)]/VPm_det;
            K=Pm*iVPm;
            
            % EKF update equations
            Pp=(eye(2)-K)*Pm;
            xp=xm+K*(y-xm);
            
            % ball norm limiter
            nrm=sqrt(Pp(:)'*Pp(:));
            Pp=Pp*min(nrm,obj.Plim)/max(nrm,realmin);
            
            % store
            obj.P = Pp;
            obj.x = xp;
            obj.vold = ve;

        end
        
        function N = getNumInputsImpl(~)
            % Specify number of System inputs
            N = 3; % Because stepImpl has one argument beyond obj
        end
        
        function N = getNumOutputsImpl(~)
            % Specify number of System outputs
            N = 1; % Because stepImpl has one output
        end
    end
end

