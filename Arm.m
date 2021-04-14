classdef Arm < handle & matlab.mixin.SetGet
   %ARM Class for simulating a two-segment viscoelastic elbow joint with rigid fixation at the shoulder.
   
   properties (Constant, Hidden, Access=public)
      dt     = 0.001; % (seconds; time-step between successive samples)
      E      = [];    % Elastic modulus [biceps, triceps]
      eta    = [];    % coefficients of viscosity [biceps, triceps]
   end
   
   properties (GetAccess=public, SetAccess=private)
      r      = [ 0.315; 0.290 ]   % Radius of [upper, lower] segments (m)
      theta  = [ 30;    15    ]   % Angle (deg) from horizontal of upper, angle at elbow
      delta  = 0.020              % How much does the elbow stick out to serve as attachment for triceps? (m)
      lambda = [0.010, 0.050; 
                0.025, 0.015]     % Connection offsets from base of segment (m) for viscoelastic muscle
      sigma     = [0, 0];            % Stress on [biceps, triceps] (N)
      d_epsilon = [0, 0];            % Change in strain on [biceps, triceps] (m)
   end
   
   properties (Hidden, Access=public)
      fig      % Figure handle
      ax       % Axes handle containing the line graphics
      lgd      % Legend handle
      upper    % Handle to line segment illustrating upper limb.
      lower    % Handle to line segment illustrating lower limb.
      biceps   % Handle to line segment illustrating biceps attachments.
      triceps  % Handle to line segment illustrating triceps attachments.
      rigid    % Handle to patch illustrating rigid attachment to body.
   end
   
   properties (Dependent, Access=public)
      tag
   end
   
   properties (Access=private)
      tag_ = '';
   end
   
   methods
      function self = Arm(name, varargin)
         %ARM Construct an instance of this class
         %
         % Examples:
         %  obj = Arm();
         %  obj = Arm('Left');  % Sets `.tag` property to 'Left'.
         %
         % Inputs
         %  name - (Optional) The 'tag' property value. 
         %  varargin - (Optional) The <'Name',Value> parameter pairs for
         %                 modifying properties in the constructor.
         if nargin < 1
            name = 'Untitled';
         end
         for iV = 1:2:numel(varargin)
            self.(varargin{iV}) = varargin{iV+1};
         end
         self.init_graphics();
         self.tag = name;
      end
      
      function delete(self)
         %DELETE Overload delete method to ensure that figure is deleted.
         if ~isempty(self.fig)
            if isvalid(self.fig)
               delete(self.fig);
            end
         end            
      end
      
      function value = get.tag(self)
         %GET.TAG Returns the "tag" (name) of the Arm object.
         value = self.tag_;
      end
      
      function set.tag(self, value)
         %SET.TAG Set the "tag" value of the Arm object.
         self.tag_ = value;
         self.fig.Name = sprintf('Arm Plot (%s)', value);
      end      
   end
   
   methods (Access=public)
      
      function configure(self,varargin)
         %CONFIGURE Updates core properties
         %
         %  self.configure('r',[values],...);
         %
         % Inputs
         %  varargin - 'Name',value parameter pairs for specific
         %                 properties:
         %              * r - radii of [upper, lower] segments (m)
         %              * theta - [degrees clockwise from horizontal for
         %                          upper segment, degrees counterclockwise
         %                          from orientation of upper segment for
         %                          lower segment]
         %              * delta - "offset" of first part (before the
         %                 "joint" attachment) for the limb segments (m)
         %              * lambda - attachment points (m) of muscles
         %         
         % This function ensures that graphics are updated to reflect
         % changes in these core parameters that affect the
         % length/relations of otherwise static simulation parameters.
         for iV = 1:2:numel(varargin)
            self.(varargin{iV}) = varargin{iV+1};
         end
         self.draw();
      end
      
      function draw(self)
         %DRAW Draw/update the graphics values.
         [x_u, y_u] = self.Upper();
         [x_l, y_l] = self.Lower();
         [x_b, y_b] = self.Biceps();
         [x_t, y_t] = self.Triceps();
         set(self.upper,'XData',x_u,'YData',y_u);
         set(self.lower,'XData',x_l,'YData',y_l);
         set(self.biceps,'XData',x_b,'YData',y_b);
         set(self.triceps,'XData',x_t,'YData',y_t);
         drawnow();
      end
      
   end
   
   methods (Hidden)
      function [x, y] = Upper(self)
         %UPPER Return the cartesian positions of upper segment
         x = zeros(1,2);
         y = zeros(1,2);
         [x(2), y(2)] = self.polar2cartesian(self.r(1), -self.theta(1));
      end
      
      function [x, y] = Lower(self)
         %LOWER Return the cartesian positions of lower segment
         x = nan(1,3);
         y = nan(1,3);
         [x(2), y(2)] = self.polar2cartesian(self.r(1), -self.theta(1));
         [a, b] = self.polar2cartesian(self.delta, self.theta(2));
         [a, b] = self.rotate(a, b, -self.theta(1));
         x(1) = x(2) - a;
         y(1) = y(2) - b;         
         [a, b] = self.polar2cartesian(self.r(2), self.theta(2));
         [a, b] = self.rotate(a, b, -self.theta(1));
         x(3) = a + x(2);
         y(3) = b + y(2);
      end
      
      function [x, y] = Biceps(self)
         %BICEPS Return the cartesian position for both biceps endpoints
         x = nan(1,2);
         y = nan(1,2);
         [x(1), y(1)] = self.polar2cartesian(self.lambda(1,1), -self.theta(1));
         [a, b] = self.polar2cartesian(self.r(1), -self.theta(1));
         [c, d] = self.polar2cartesian(self.lambda(2,1), self.theta(2));
         [c, d] = self.rotate(c, d, -self.theta(1));
         x(2) = a + c;
         y(2) = b + d;
      end
      
      function [x, y] = Triceps(self)
         %TRICEPS Return the cartesian position for both triceps endpoints
         x = nan(1,2);
         y = nan(1,2);
         [x(1), y(1)] = self.polar2cartesian(self.lambda(1,2), -self.theta(1));
         [a, b] = self.polar2cartesian(self.r(1), -self.theta(1));
         [c, d] = self.polar2cartesian(self.lambda(2,2), self.theta(2));
         [c, d] = self.rotate(c, d, -self.theta(1));
         x(2) = a - c;
         y(2) = b - d;
      end
         
   end
   
   methods (Access=private)
      function init_graphics(self)
         %INIT_GRAPHICS Initialize the figure graphics at onset.
         f = figure('Name', 'Arm Plot', 'Color', 'k', ...
            'Units', 'Normalized', 'Position', [0.2 0.2 0.4 0.4], ...
            'Visible', true);
         h = axes(f, 'NextPlot', 'add', 'Color', 'k', ...
            'LineWidth', 1.5, 'XColor','w', 'YColor', 'w',...
            'FontName','Arial','FontSize',13,...
            'XLim',[-0.1, 0.4], 'XTick',[0 0.2 0.4], ...
            'YLim',[-0.25 0.25], 'YTick',[-0.1 0 0.1]);
         xlabel(h,'AP (m)','FontName','Arial','Color','w','FontSize',14);
         ylabel(h,'ML (m)','FontName','Arial','Color','w','FontSize',14);
         title(h,'Ball-and-Stick 2-Segment Arm', ...
            'FontName','Arial','Color','w','FontSize',16,'FontWeight','bold');
         self.rigid = fill(h, ...
            [-0.2, 0, 0, -0.2], ... % XData
            [-1 -1 1 1], ... % YData
            validatecolor('#333'), ... % Color
            'FaceColor',validatecolor('#333'),... % Color
            'EdgeColor',validatecolor('#888'),...
            'DisplayName','Body');
         [x_u, y_u] = self.Upper();
         [x_l, y_l] = self.Lower();
         [x_b, y_b] = self.Biceps();
         [x_t, y_t] = self.Triceps();
         
         self.upper = line(h, x_u, y_u,...
            'DisplayName','Upper Arm', ...
            'LineWidth', 5,...
            'Color','#d5dbcc',...
            'MarkerSize', 10, ...
            'Marker', 's', ...
            'MarkerFaceColor', '#555', ...
            'MarkerEdgeColor', 'k', ...
            'MarkerIndices', 1);
         
         self.lower = line(h, x_l, y_l,...
            'DisplayName','Lower Arm', ...
            'LineWidth', 5,...
            'Color','#d5dbcc',...
            'MarkerSize', 8, ...
            'Marker', 'o', ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', 'r', ...
            'MarkerIndices', 2);
         
         self.biceps = line(h, x_b, y_b,...
            'DisplayName','Biceps', ...
            'LineWidth', 3.5,...
            'LineStyle',':', ...
            'Color','#d2340f',...
            'MarkerSize', 8, ...
            'Marker', 's', ...
            'MarkerFaceColor', '#c0d1a9', ...
            'MarkerEdgeColor', '#666');
         
         self.triceps = line(h, x_t, y_t,...
            'DisplayName','Triceps', ...
            'LineWidth', 3.5,...
            'LineStyle',':', ...
            'Color',' #d2870f',...
            'MarkerSize', 8, ...
            'Marker', 's', ...
            'MarkerFaceColor', '#c0d1a9', ...
            'MarkerEdgeColor', '#666');
         
         self.lgd = legend(h, ...
            'Location', 'Northeast', ...
            'FontName','Arial', ...
            'Color','none',...
            'EdgeColor','none',...
            'TextColor','w',...
            'FontSize', 12);
         
         self.fig = f;
         self.ax = h;
      end
   end
   
   methods (Static)
      function [x, y] = polar2cartesian(r, theta)
         %POLAR2CARTESIAN Convert polar coordinates to cartesian.
         x = r * cos(deg2rad(theta)); 
         y = r * sin(deg2rad(theta));
      end
      
      function [x, y] = rotate(x, y, theta)
         %ROTATE Rotate cartesian position by theta degrees, clockwise is
         %        positive direction.
         theta = deg2rad(theta);
         pos = [cos(theta), sin(theta);
                -sin(theta), cos(theta)] * [x; y];
         x = pos(1);
         y = pos(2);             
      end
      
      function test(a,b,varargin)
         disp(a);
         for iV = 1:numel(varargin)
            disp(varargin{iV});
         end
         disp(b);
      end
   end
end

