function plot_error_convergence(h, err)
    % 验证输入数据
    if nargin ~= 2
        error('需要两个输入参数: h 和 err');
    end
    if ~isvector(h) || ~isvector(err) || length(h) ~= length(err)
        error('h 和 err 必须是相同长度的向量');
    end
    if any(h <= 0) || any(err <= 0)
        error('h 和 err 必须包含正数');
    end
    
    % 计算对数坐标
    x = -log(h(:));  % 确保列向量
    y = -log(err(:));
    
    % 最小二乘线性拟合
    A = [x, ones(size(x))];
    coeffs = A \ y;  % 求解线性方程组
    slope = coeffs(1);
    intercept = coeffs(2);
    
    % 生成拟合曲线数据
    x_fit = linspace(min(x), max(x), 100);
    y_fit = slope * x_fit + intercept;
    
    % 创建图形窗口
    fig = figure('Name', '误差收敛图 - 可拖动公式', 'NumberTitle', 'off', ...
                 'Position', [100, 100, 800, 600], 'Color', 'w');
    
    % 绘制原始数据点和拟合曲线
    ax = axes('Parent', fig, 'FontSize', 11);
    scatter(ax, x, y, 70, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8], ...
            'MarkerEdgeColor', 'k', 'DisplayName', '原始数据');
    hold(ax, 'on');
    plot(ax, x_fit, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', '最小二乘拟合');
    grid(ax, 'on');
    
    % 设置坐标轴标签和标题
    xlabel(ax, '$-\log(h)$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(ax, '$-\log(err)$', 'Interpreter', 'latex', 'FontSize', 12);
    title(ax, '误差收敛图', 'FontSize', 14);
    
    % 计算数据范围
    x_range = range(x);
    y_range = range(y);
    
    % 初始公式位置（曲线中点上方）
    mid_idx = round(0.5 * length(x_fit));
    x_mid = x_fit(mid_idx);
    y_mid = y_fit(mid_idx);
    x_offset = 0.05 * x_range;
    y_offset = 0.1 * y_range;
    
    % 创建可拖动的公式文本框
    eqn = sprintf('$y = %.4f x + %.4f$', slope, intercept);
    formula_text = text(ax, x_mid + x_offset, y_mid + y_offset, eqn, ...
        'Interpreter', 'latex', 'FontSize', 12, ...
        'BackgroundColor', [1, 1, 1, 0.8], ...  % 半透明白色背景
        'EdgeColor', 'k', 'Margin', 3, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'ButtonDownFcn', @startDrag, ...
        'Tag', 'DraggableText');
    
    % 添加图例
    legend(ax, 'Location', 'northwest');
    
    % 添加说明文本
    uicontrol(fig, 'Style', 'text', 'String', '提示：点击并拖动公式可自由移动位置', ...
              'Position', [10, 10, 300, 20], 'BackgroundColor', [1, 1, 1, 0.7], ...
              'FontSize', 10, 'HorizontalAlignment', 'left');
    
    % 调整坐标轴范围确保所有元素可见
    xlim(ax, [min(x)-0.1*x_range, max(x)+0.2*x_range]);
    ylim(ax, [min(y)-0.1*y_range, max(y)+0.2*y_range]);
    
    % 存储拖动状态
    setappdata(fig, 'isDragging', false);
    setappdata(fig, 'dragOffset', [0, 0]);
    
    % 设置窗口回调函数
    set(fig, 'WindowButtonMotionFcn', @dragging, ...
             'WindowButtonUpFcn', @stopDrag);
    
    % 拖动开始函数
    function startDrag(src, ~)
        setappdata(fig, 'isDragging', true);
        
        % 获取当前点位置
        currentPoint = get(ax, 'CurrentPoint');
        textPos = get(src, 'Position');
        
        % 计算偏移量
        offset = [textPos(1) - currentPoint(1,1), textPos(2) - currentPoint(1,2)];
        setappdata(fig, 'dragOffset', offset);
        
        % 高亮显示被拖动的文本
        set(src, 'BackgroundColor', [0.9, 0.95, 1, 0.9], 'EdgeColor', 'b', 'LineWidth', 1.5);
    end

    % 拖动中函数
    function dragging(~, ~)
        if getappdata(fig, 'isDragging')
            currentPoint = get(ax, 'CurrentPoint');
            offset = getappdata(fig, 'dragOffset');
            
            % 更新文本位置
            newPos = [currentPoint(1,1) + offset(1), currentPoint(1,2) + offset(2)];
            set(formula_text, 'Position', [newPos(1), newPos(2), 0]);
            
            % 强制刷新图形
            drawnow;
        end
    end

    % 拖动结束函数
    function stopDrag(~, ~)
        if getappdata(fig, 'isDragging')
            setappdata(fig, 'isDragging', false);
            
            % 恢复文本样式
            set(formula_text, 'BackgroundColor', [1, 1, 1, 0.8], ...
                'EdgeColor', 'k', 'LineWidth', 0.5);
        end
    end
end