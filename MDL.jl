using Markdown
using Makie


const Weights = [1/8, 1/4, 1/2]


@doc md"""
    用于存储不同离子的数据
""" ->
function _IonData(choice::Int64)
    Data::Array{Tuple{Float64, Float64, Array{Float64}}} = [
        #=
            对于每一种离子按照顺序记录以下值
                自身电荷、第一种层离子（该层离子均等同），第二种层离子
            离子种类按照
                C_60、K(1)、K(2)
            的顺序记录
        =#
        (-3, 1, [1, -3]),
        (1, 1, [-3, 1]),
        (1, -1, [1, 1]) # 这里是将C_60和K(1)对K(2)的电磁作用平均（反正对称）
    ]

    return Data[choice]
end


@doc md"""
    最近的离子与参考点的距离为 $4$
    此为晶格常数，将不同的离子分层计算并输出结果成图像
""" ->
function _Calculate(n::Int64, choice::Int64)
    e::Float64, layer1::Float64, Layer2::Array{Float64} = _IonData(choice)
    name::String = ["C60", "K(1)", "K(2)"][choice]

    Points::Array{Float64} = []
    Edges::Array{Float64} = [0]
    Faces::Array{Float64} = [0]

    # 处理顶角离子
    for i::Int64 = 1:n
        point::Float64 = 8/sqrt(3i^2)
        if mod(i, 4) == 0
            # 顶角回到参考点离子，完成一个正常的周期
            point *= e
        else
            #=
                如果为奇数则为第一种层离子
                如果为偶数则为第二种层中固定的顶角离子
            =#
            isodd(i) && (point *= layer1; true) || (point *= Layer2[1]; true)
        end
        push!(Points, point)
    end
    # 处理棱离子
    #=
        从第二层开始出现棱离子
        如果层数为奇数则为第一种层离子
        如果层数为偶数则为第二种层离子
    =#
    for i::Int64=2:n
        edge::Float64 = 0
        if iseven(i)
            for j::Int64 = 2-i:2:i-2
                temp::Float64 = 12/sqrt(2i^2+j^2)
                isodd(j ÷ 2) && (temp *= Layer2[1]; true) || (temp *= Layer2[2]; true)
                edge += temp
            end
        else
            for j::Int64 = 1:2:i-2
                temp::Float64 = 24/sqrt(2i^2 + j^2)
                temp *= layer1
                edge += temp
            end
        end
        push!(Edges, edge)
    end
    # 处理面离子
    #=
        从第二层开始出现面离子
        如果层数为奇数则为第一种层离子
        如果层数为偶数则为第二种层离子
    =#
    for i::Int64 = 2:n
        face::Float64 = 0
        if iseven(i)
            for j::Int64 = 2-i:2:i-2
                for k::Int64 = 2-i:2:i-2
                    temp::Float64 = 6/sqrt(i^2 + j^2 + k^2)
                    if mod(i,4) == 0
                        isodd(j ÷ 2) != isodd(k÷2) && (temp *= Layer2[1]; true) || (temp*=Layer2[2]; true)
                    else
                        isodd(j ÷ 2) != isodd(k÷2) && (temp *= Layer2[2]; true) || (temp*=Layer2[1]; true)
                    end
                    face += temp
                end
            end
        else
            for j::Int64 = 1:2:i-2
                for k::Int64 = 1:2:i-2
                    temp::Float64 = 24/sqrt(i^2 + j^2 + k^2)
                    temp *= layer1
                    face += temp
                end
            end
        end
        push!(Faces, face)
    end

    Ions::Array{Array{Float64}} = [Points, Edges, Faces]

    Α::Array{Float64} = []
    α::Float64 = 0
    β::Float64 = 0
    for i=1:n
        Temp = Ions .|> (X -> X[i])
        α = β + sum(Weights .* Temp)
        β += Temp |> sum
        push!(Α,α)
    end
    Α .*= -e * 4

    FixΑ::Array{Float64} = @. (Α[1:end-3] + Α[2:end-2] + Α[3:end-1] + Α[4:end])/4

    FixRange = range(1, n, length = length(FixΑ))
    Range = 1:n
    scene = lines(FixRange, FixΑ, color = :blue)
    scatter!(FixRange, FixΑ, color = :blue, markersize = 1)
    scene = title("处理后$(name)的马德隆常数")
    Makie.save("Fixed-$(name).png", scene)
    scene = scatter(Range, Α,color = :red, markersize = 1)
    scene = title("未经处理$(name)的马德隆常数")
    Makie.save("Raw-$(name).png", scene)

    return FixΑ[end]
end


@doc md"""
    用于确定晶格个数、测算时间消耗以及得出最后结果
""" ->
function _Main()
    Names::Array{String} = ["C60", "K(1)", "K(2)"]

    println("请输入所需晶格个数：")
    n = parse(Int64, readline(stdin))
    println()
    
    for i = 1:3
        name::String = Names[i]
        println("开始计算$(name)的马德隆常数")
        println("*"^10)

        α::Float64 = 0
        @time α=_Calculate(n, i)
        println("最终马德隆常数为 $(α)")
        println("*"^10)
        println()
    end
end
