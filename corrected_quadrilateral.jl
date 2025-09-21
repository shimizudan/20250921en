using LinearAlgebra

function max_area_convex_quadrilateral_corrected(a, b, c, d)
    """
    修正版：4辺の長さから面積最大の凸四角形を求める
    理論上限の計算を修正
    """
    
    # まず四角形が構成可能かチェック（三角不等式）
    sides = [a, b, c, d]
    for i in 1:4
        sum_others = sum(sides) - sides[i]
        if sides[i] >= sum_others
            println("エラー: 三角不等式を満たしません")
            return nothing
        end
    end
    
    # ピトー偏差の計算
    pitot_deviation = abs((a + c) - (b + d))
    println("ピトー偏差: $(round(pitot_deviation, digits=4))")
    
    # 理論上限の正しい計算
    if pitot_deviation < 1e-6
        # ピトー条件を満たす場合：ブラーマグプタ公式が適用可能
        s = (a + b + c + d) / 2
        discriminant = (s - a) * (s - b) * (s - c) * (s - d)
        if discriminant > 0
            theoretical_max = sqrt(discriminant)
            println("理論上限面積（ブラーマグプタ）: $(round(theoretical_max, digits=6))")
        else
            println("エラー: ブラーマグプタの判別式が負です")
            return nothing
        end
    else
        # ピトー条件を満たさない場合：理論上限は不明
        println("注意: ピトー条件を満たさないため、ブラーマグプタ公式は適用不可")
        
        # 代替的な上限の推定（対角線を使った面積の上限）
        # 最大可能対角線長を使った推定
        max_diag1 = a + c  # 最大対角線1
        max_diag2 = b + d  # 最大対角線2
        
        # 対角線が直交する場合の面積上限
        estimated_upper_bound = max_diag1 * max_diag2 / 2
        println("推定上限面積（対角線直交時）: $(round(estimated_upper_bound, digits=6))")
        theoretical_max = estimated_upper_bound
    end
    
    best_area = 0.0
    best_result = nothing
    
    println("グリッドサーチによる最適化を実行中...")
    
    # 内角ベースのグリッドサーチ
    angle_steps = 15
    count = 0
    total = angle_steps^3
    
    for i1 in 1:angle_steps
        for i2 in 1:angle_steps
            for i3 in 1:angle_steps
                count += 1
                if count % 100 == 0
                    progress = round(count/total*100, digits=1)
                    print("進捗: $(progress)%\r")
                end
                
                # 内角を生成（20度から160度の範囲）
                α = π/9 + (8*π/9) * (i1 - 1) / (angle_steps - 1)  # 20°-160°
                β = π/9 + (8*π/9) * (i2 - 1) / (angle_steps - 1)
                γ = π/9 + (8*π/9) * (i3 - 1) / (angle_steps - 1)
                δ = 2π - α - β - γ
                
                # 4番目の角度が有効範囲内かチェック
                if δ < π/9 || δ > 8*π/9
                    continue
                end
                
                # この角度で四角形を構築
                vertices = construct_quadrilateral_robust(a, b, c, d, α, β, γ, δ)
                
                if vertices !== nothing
                    if is_convex_and_valid_sides(vertices, [a, b, c, d])
                        area = polygon_area(vertices)
                        if area > best_area
                            best_area = area
                            best_result = (vertices, area, theoretical_max, [α, β, γ, δ])
                        end
                    end
                end
            end
        end
    end
    
    println("\n最適化完了")
    
    if best_result === nothing
        println("有効な解が見つかりませんでした")
        return nothing
    end
    
    return best_result
end

function construct_quadrilateral_robust(a, b, c, d, α, β, γ, δ)
    """
    より確実な四角形構築アルゴリズム
    """
    
    # P1を原点、P2を正のx軸上に固定
    P1 = [0.0, 0.0]
    P2 = [a, 0.0]
    
    # P3の位置を内角βから計算
    # P2での内角がβの場合、P2からP3への角度は π - β
    angle_P2_to_P3 = π - β
    P3 = [P2[1] + b * cos(angle_P2_to_P3), 
          P2[2] + b * sin(angle_P2_to_P3)]
    
    # P4の位置を2つの制約から求める
    # 制約1: |P1 - P4| = d
    # 制約2: |P3 - P4| = c
    
    # 2つの円の交点を求める
    P4_candidates = circle_intersection_robust(P1, d, P3, c)
    
    for P4 in P4_candidates
        vertices = [P1, P2, P3, P4]
        
        # 辺の長さをチェック
        actual_sides = [
            norm(P2 - P1),  # a
            norm(P3 - P2),  # b
            norm(P4 - P3),  # c
            norm(P1 - P4)   # d
        ]
        
        # 辺の長さの誤差をチェック
        side_errors = abs.([a, b, c, d] - actual_sides)
        
        if all(error -> error < 0.01, side_errors)
            # 内角もチェック
            if verify_angles(vertices, [α, β, γ, δ])
                return vertices
            end
        end
    end
    
    return nothing
end

function verify_angles(vertices, target_angles)
    """
    四角形の内角が目標角度と一致するかチェック
    """
    n = length(vertices)
    actual_angles = []
    
    for i in 1:n
        p1 = vertices[i]
        p2 = vertices[i % n + 1]
        p3 = vertices[(i + 1) % n + 1]
        
        # ベクトル
        v1 = p1 - p2
        v2 = p3 - p2
        
        # 内角の計算
        cos_angle = dot(v1, v2) / (norm(v1) * norm(v2))
        cos_angle = clamp(cos_angle, -1.0, 1.0)
        angle = acos(cos_angle)
        
        push!(actual_angles, angle)
    end
    
    # 角度の誤差をチェック
    angle_errors = abs.(actual_angles - target_angles)
    
    return all(error -> error < 0.1, angle_errors)  # 約6度の許容誤差
end

function circle_intersection_robust(center1, radius1, center2, radius2)
    """
    2つの円の交点を求める（堅牢版）
    """
    d = norm(center2 - center1)
    
    # 円が交わらない場合
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
    
    if h_squared < -1e-10
        return []
    end
    
    h = sqrt(max(0, h_squared))
    
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

function is_convex_and_valid_sides(vertices, target_sides)
    """
    凸性と辺の長さの両方をチェック
    """
    
    # 凸性チェック
    if !is_convex_polygon_robust(vertices)
        return false
    end
    
    # 辺の長さチェック
    actual_sides = []
    n = length(vertices)
    
    for i in 1:n
        j = i % n + 1
        side_length = norm(vertices[j] - vertices[i])
        push!(actual_sides, side_length)
    end
    
    side_errors = abs.(actual_sides - target_sides)
    
    return all(error -> error < 0.05, side_errors)
end

function is_convex_polygon_robust(vertices)
    """
    堅牢な凸性チェック
    """
    n = length(vertices)
    if n < 3
        return false
    end
    
    cross_products = []
    
    for i in 1:n
        p1 = vertices[i]
        p2 = vertices[i % n + 1]
        p3 = vertices[(i + 1) % n + 1]
        
        v1 = p2 - p1
        v2 = p3 - p2
        
        cross_product = v1[1] * v2[2] - v1[2] * v2[1]
        
        if abs(cross_product) > 1e-10  # 数値誤差を考慮
            push!(cross_products, sign(cross_product))
        end
    end
    
    # すべて同じ符号かチェック
    return length(unique(cross_products)) <= 1
end

function polygon_area(vertices)
    """
    Shoelace公式による面積計算
    """
    n = length(vertices)
    area = 0.0
    
    for i in 1:n
        j = i % n + 1
        area += vertices[i][1] * vertices[j][2]
        area -= vertices[j][1] * vertices[i][2]
    end
    
    return abs(area) / 2
end

# Bretschneider公式による正確な計算
function bretschneider_area(a, b, c, d, p, q, φ)
    """
    Bretschneider公式による面積計算
    p, q: 対角線の長さ
    φ: 対角線の交差角
    """
    
    # Bretschneider公式: K = √[(s-a)(s-b)(s-c)(s-d) - abcd·cos²((α+γ)/2)]
    # ここで α+γ は対向する角の和
    
    s = (a + b + c + d) / 2
    
    # 対角線の関係から対向角の和を計算
    # cos(α+γ) を対角線の長さから求める
    cos_alpha_plus_gamma = (p^2 + q^2 - a^2 - c^2) / (2 * sqrt((a^2 + b^2) * (c^2 + d^2)))
    
    area_squared = (s-a)*(s-b)*(s-c)*(s-d) - a*b*c*d*cos_alpha_plus_gamma^2
    
    if area_squared < 0
        return 0
    end
    
    return sqrt(area_squared)
end

# デモ関数
function demo_corrected()
    println("=== 修正版：四角形面積最大化（理論上限修正版） ===")
    
    test_cases = [
        (3.0, 4.0, 5.0, 6.0),
        (4.0, 5.0, 4.0, 5.0),  # ピトー条件を満たす
        (3.0, 4.0, 4.0, 3.0),  # ピトー条件を満たす
        (5.0, 5.0, 5.0, 5.0),  # 正方形
        (2.0, 3.0, 4.0, 5.0)
    ]
    
    for (a, b, c, d) in test_cases
        println("\n" * "="^60)
        println("辺の長さ: a=$a, b=$b, c=$c, d=$d")
        
        result = max_area_convex_quadrilateral_corrected(a, b, c, d)
        
        if result !== nothing
            vertices, area, theoretical_max, angles = result
            
            println("\n最終結果:")
            println("最大面積: $(round(area, digits=6))")
            println("推定上限: $(round(theoretical_max, digits=6))")
            
            # 効率の計算（ピトー条件を満たす場合のみ意味がある）
            pitot_deviation = abs((a + c) - (b + d))
            if pitot_deviation < 1e-6
                efficiency = area / theoretical_max * 100
                println("効率（ブラーマグプタに対する）: $(round(efficiency, digits=2))%")
            else
                ratio = area / theoretical_max * 100
                println("対推定上限比: $(round(ratio, digits=2))%")
            end
            
            println("\n最適内角（度）:")
            for (i, angle) in enumerate(angles)
                println("  角$i: $(round(rad2deg(angle), digits=2))°")
            end
            
            println("\n頂点座標:")
            for (i, v) in enumerate(vertices)
                println("  P$i: ($(round(v[1], digits=4)), $(round(v[2], digits=4)))")
            end
            
            # 検証
            actual_sides = []
            for i in 1:4
                j = i % 4 + 1
                push!(actual_sides, norm(vertices[j] - vertices[i]))
            end
            println("\n検証:")
            println("実際の辺の長さ: $(round.(actual_sides, digits=4))")
            println("目標の辺の長さ: [$a, $b, $c, $d]")
            println("凸四角形: $(is_convex_polygon_robust(vertices))")
            
            # ピトー偏差の再確認
            println("ピトー偏差 |a+c-b-d|: $(round(abs((a + c) - (b + d)), digits=4))")
        else
            println("解が見つかりませんでした")
        end
    end
end

# 実行
demo_corrected()