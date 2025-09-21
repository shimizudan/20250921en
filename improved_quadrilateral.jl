using LinearAlgebra

function max_area_convex_quadrilateral_improved(a, b, c, d)
    """
    改良版：4辺の長さから面積最大の凸四角形を求める
    より確実な最適化アルゴリズムを使用
    """
    
    # 理論上限（ブラーマグプタ公式）
    s = (a + b + c + d) / 2
    discriminant = (s - a) * (s - b) * (s - c) * (s - d)
    
    if discriminant <= 0
        println("エラー: この辺の長さでは四角形を構成できません")
        return nothing
    end
    
    theoretical_max = sqrt(discriminant)
    println("理論上限面積: $(round(theoretical_max, digits=6))")
    
    # 最適化変数：対角線の長さと角度
    best_area = 0.0
    best_result = nothing
    
    # 対角線ベースのアプローチ
    println("対角線ベース最適化を実行中...")
    
    # 対角線の長さの範囲を設定
    min_diag = max(abs(a-c), abs(b-d)) + 0.1
    max_diag = (a+c+b+d) - 0.1
    
    diag_steps = 30
    angle_steps = 20
    
    count = 0
    total = diag_steps * diag_steps * angle_steps
    
    for i1 in 1:diag_steps
        for i2 in 1:diag_steps
            for i3 in 1:angle_steps
                count += 1
                if count % 500 == 0
                    progress = round(count/total*100, digits=1)
                    print("進捗: $(progress)%\r")
                end
                
                # 対角線の長さ
                p = min_diag + (max_diag - min_diag) * (i1 - 1) / (diag_steps - 1)
                q = min_diag + (max_diag - min_diag) * (i2 - 1) / (diag_steps - 1)
                
                # 対角線の交点における角度
                θ = π/6 + (5π/6 - π/6) * (i3 - 1) / (angle_steps - 1)
                
                # この対角線で四角形を構築
                vertices = construct_quadrilateral_from_diagonals(a, b, c, d, p, q, θ)
                
                if vertices !== nothing
                    if is_convex_and_valid(vertices, [a, b, c, d])
                        area = polygon_area(vertices)
                        if area > best_area
                            best_area = area
                            best_result = (vertices, area, theoretical_max, [p, q, θ])
                        end
                    end
                end
            end
        end
    end
    
    println("\n最適化完了")
    
    if best_result === nothing
        println("有効な解が見つかりませんでした")
        # フォールバック：シンプレックス法
        return simplex_optimization(a, b, c, d, theoretical_max)
    end
    
    return best_result
end

function construct_quadrilateral_from_diagonals(a, b, c, d, p, q, θ)
    """
    対角線の長さと交差角から四角形を構築
    p, q: 対角線の長さ
    θ: 対角線の交差角
    """
    
    # 対角線の交点を原点とする
    O = [0.0, 0.0]
    
    # 対角線の分割比を変数として扱う
    for t1 in 0.2:0.1:0.8  # P1からの分割比
        for t2 in 0.2:0.1:0.8  # P2からの分割比
            
            # 頂点位置の計算
            P1 = [-p * t1, 0.0]
            P3 = [p * (1 - t1), 0.0]
            P2 = [q * t2 * cos(θ), q * t2 * sin(θ)]
            P4 = [-q * (1 - t2) * cos(θ), -q * (1 - t2) * sin(θ)]
            
            vertices = [P1, P2, P3, P4]
            
            # 辺の長さチェック
            sides = calculate_side_lengths(vertices)
            target_sides = [a, b, c, d]
            
            errors = abs.(sides - target_sides)
            
            if all(e -> e < 0.05, errors)  # 許容誤差内
                return vertices
            end
        end
    end
    
    return nothing
end

function calculate_side_lengths(vertices)
    """頂点から辺の長さを計算"""
    n = length(vertices)
    sides = []
    
    for i in 1:n
        j = i % n + 1
        side_length = norm(vertices[j] - vertices[i])
        push!(sides, side_length)
    end
    
    return sides
end

function is_convex_and_valid(vertices, target_sides)
    """凸性と辺の長さをチェック"""
    
    # 凸性チェック
    if !is_convex_polygon(vertices)
        return false
    end
    
    # 辺の長さチェック
    actual_sides = calculate_side_lengths(vertices)
    errors = abs.(actual_sides - target_sides)
    
    return all(e -> e < 0.1, errors)
end

function is_convex_polygon(vertices)
    """改良された凸性チェック"""
    n = length(vertices)
    if n < 3
        return false
    end
    
    signs = []
    
    for i in 1:n
        p1 = vertices[i]
        p2 = vertices[i % n + 1]
        p3 = vertices[(i + 1) % n + 1]
        
        # 外積の計算
        v1 = p2 - p1
        v2 = p3 - p2
        cross = v1[1] * v2[2] - v1[2] * v2[1]
        
        if abs(cross) > 1e-10  # 数値誤差を考慮
            push!(signs, sign(cross))
        end
    end
    
    # すべて同じ符号かチェック
    return length(unique(signs)) <= 1
end

function polygon_area(vertices)
    """Shoelace公式による面積計算"""
    n = length(vertices)
    area = 0.0
    
    for i in 1:n
        j = i % n + 1
        area += vertices[i][1] * vertices[j][2]
        area -= vertices[j][1] * vertices[i][2]
    end
    
    return abs(area) / 2
end

function simplex_optimization(a, b, c, d, theoretical_max)
    """
    フォールバック：シンプレックス法による最適化
    """
    println("シンプレックス法による最適化を開始...")
    
    best_area = 0.0
    best_result = nothing
    
    # ランダムサーチ + 局所最適化
    for trial in 1:1000
        # ランダムな初期値
        p = rand() * (a + c) + abs(a - c)/2
        q = rand() * (b + d) + abs(b - d)/2
        θ = π/6 + rand() * (2π/3)
        
        # 局所最適化
        for iter in 1:10
            for δ in [-0.1, 0.1]
                for var in 1:3
                    if var == 1
                        test_p = p + δ
                        test_q = q
                        test_θ = θ
                    elseif var == 2
                        test_p = p
                        test_q = q + δ
                        test_θ = θ
                    else
                        test_p = p
                        test_q = q
                        test_θ = θ + δ/10
                    end
                    
                    vertices = construct_quadrilateral_from_diagonals(a, b, c, d, test_p, test_q, test_θ)
                    
                    if vertices !== nothing && is_convex_and_valid(vertices, [a, b, c, d])
                        area = polygon_area(vertices)
                        if area > best_area
                            best_area = area
                            best_result = (vertices, area, theoretical_max, [test_p, test_q, test_θ])
                            p, q, θ = test_p, test_q, test_θ
                        end
                    end
                end
            end
        end
    end
    
    return best_result
end

# より効率的なバージョン：制約付き最適化
function max_area_constrained_optimization(a, b, c, d)
    """
    制約付き最適化による解法
    """
    
    s = (a + b + c + d) / 2
    discriminant = (s - a) * (s - b) * (s - c) * (s - d)
    
    if discriminant <= 0
        println("エラー: この辺の長さでは四角形を構成できません")
        return nothing
    end
    
    theoretical_max = sqrt(discriminant)
    println("理論上限面積: $(round(theoretical_max, digits=6))")
    
    # Bretschneider公式を使用した直接的アプローチ
    best_area = 0.0
    best_result = nothing
    
    println("Bretschneider公式による最適化...")
    
    # 対角線の角度を変数とする
    for i in 1:100
        φ = π * i / 100  # 対角線の角度
        
        # Bretschneider公式から面積を計算
        cos_phi = cos(φ)
        
        # (a²+c²-b²-d²)cos(φ) の係数
        coeff = (a^2 + c^2 - b^2 - d^2)
        
        area_squared = (s-a)*(s-b)*(s-c)*(s-d) - a*b*c*d*cos_phi^2
        
        if area_squared > 0
            area = sqrt(area_squared)
            
            if area > best_area
                best_area = area
                
                # この角度に対応する四角形を構築
                vertices = construct_from_bretschneider(a, b, c, d, φ)
                
                if vertices !== nothing
                    best_result = (vertices, area, theoretical_max, [φ])
                end
            end
        end
    end
    
    return best_result
end

function construct_from_bretschneider(a, b, c, d, φ)
    """
    Bretschneider公式の角度から四角形を構築
    """
    
    # コサインの法則を使用して対角線の長さを計算
    p_squared = a^2 + c^2 - 2*a*c*cos(π - φ)
    q_squared = b^2 + d^2 - 2*b*d*cos(φ)
    
    if p_squared <= 0 || q_squared <= 0
        return nothing
    end
    
    p = sqrt(p_squared)
    q = sqrt(q_squared)
    
    # 四角形を構築
    # P1を原点、P2をx軸上に配置
    P1 = [0.0, 0.0]
    P2 = [a, 0.0]
    
    # P3の位置をコサインの法則から計算
    cos_angle_p1 = (a^2 + p^2 - c^2) / (2*a*p)
    
    if abs(cos_angle_p1) > 1
        return nothing
    end
    
    angle_p1 = acos(cos_angle_p1)
    P3 = [p * cos(angle_p1), p * sin(angle_p1)]
    
    # P4の位置を制約から計算
    # |P2-P4| = b, |P3-P4| = d
    P4_candidates = circle_intersection(P2, b, P3, d)
    
    for P4 in P4_candidates
        vertices = [P1, P2, P4, P3]  # 順序を調整
        
        if is_convex_polygon(vertices)
            actual_sides = calculate_side_lengths(vertices)
            target_sides = [a, b, c, d]
            errors = abs.(actual_sides - target_sides)
            
            if all(e -> e < 0.05, errors)
                return vertices
            end
        end
    end
    
    return nothing
end

function circle_intersection(center1, radius1, center2, radius2)
    """2つの円の交点を求める（改良版）"""
    d = norm(center2 - center1)
    
    # 交わらない場合
    if d > radius1 + radius2 + 1e-10 || d < abs(radius1 - radius2) - 1e-10
        return []
    end
    
    # 同心円の場合
    if d < 1e-10
        return []
    end
    
    # 交点計算
    a = (radius1^2 - radius2^2 + d^2) / (2 * d)
    h_squared = radius1^2 - a^2
    
    if h_squared < 0
        return []
    end
    
    h = sqrt(h_squared)
    
    # 中点
    px = center1[1] + a * (center2[1] - center1[1]) / d
    py = center1[2] + a * (center2[2] - center1[2]) / d
    
    if h < 1e-10  # 接触の場合
        return [[px, py]]
    end
    
    # 2つの交点
    intersection1 = [px + h * (center2[2] - center1[2]) / d,
                     py - h * (center2[1] - center1[1]) / d]
    
    intersection2 = [px - h * (center2[2] - center1[2]) / d,
                     py + h * (center2[1] - center1[1]) / d]
    
    return [intersection1, intersection2]
end

# デモ関数
function demo_improved()
    println("=== 改良版：交差しない四角形の面積最大化 ===")
    
    test_cases = [
        (3.0, 4.0, 5.0, 6.0),
        (4.0, 5.0, 4.0, 5.0),
        (3.0, 4.0, 4.0, 3.0),
        (5.0, 5.0, 5.0, 5.0),
        (2.0, 3.0, 4.0, 5.0)
    ]
    
    for (a, b, c, d) in test_cases
        println("\n" * "="^60)
        println("辺の長さ: a=$a, b=$b, c=$c, d=$d")
        
        # 複数の手法を試す
        println("\n方法1: 改良版対角線ベース最適化")
        result1 = max_area_convex_quadrilateral_improved(a, b, c, d)
        
        println("\n方法2: Bretschneider公式ベース最適化")
        result2 = max_area_constrained_optimization(a, b, c, d)
        
        # 最良の結果を選択
        best_result = nothing
        if result1 !== nothing && result2 !== nothing
            best_result = result1[2] > result2[2] ? result1 : result2
        elseif result1 !== nothing
            best_result = result1
        elseif result2 !== nothing
            best_result = result2
        end
        
        if best_result !== nothing
            vertices, area, theoretical_max, params = best_result
            
            println("\n最終結果:")
            println("最大面積: $(round(area, digits=6))")
            println("理論上限: $(round(theoretical_max, digits=6))")
            println("効率: $(round(area/theoretical_max*100, digits=2))%")
            
            println("\n頂点座標:")
            for (i, v) in enumerate(vertices)
                println("  P$i: ($(round(v[1], digits=4)), $(round(v[2], digits=4)))")
            end
            
            # 検証
            actual_sides = calculate_side_lengths(vertices)
            println("\n検証:")
            println("実際の辺の長さ: $(round.(actual_sides, digits=4))")
            println("目標の辺の長さ: [$a, $b, $c, $d]")
            println("凸四角形: $(is_convex_polygon(vertices))")
            
            # ピトー偏差
            pitot_deviation = abs((a + c) - (b + d))
            println("ピトー偏差: $(round(pitot_deviation, digits=4))")
        else
            println("解が見つかりませんでした")
        end
    end
end

# 実行
demo_improved()