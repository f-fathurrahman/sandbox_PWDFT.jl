using Printf

function parse_namelist(filename::AbstractString)
    lines = String[]
    open(filename) do f
        for line in eachline(f)
            if occursin("!", line)
                line = split(line, "!")[1]
            end
            line = strip(line)
            if !isempty(line)
                push!(lines, line)
            end
        end
    end

    groups = Dict{String, Dict{String, Any}}()
    i = 1
    n = length(lines)

    while i <= n
        line = lines[i]
        if startswith(line, "&")
            group_name = strip(line[2:end])
            group_lines = String[]
            i += 1
            while i <= n
                l = lines[i]
                if occursin("/", l)
                    parts = split(l, "/")
                    if !isempty(parts[1])
                        push!(group_lines, parts[1])
                    end
                    break
                else
                    push!(group_lines, l)
                    i += 1
                end
            end
            groups[group_name] = parse_group(group_lines)
        end
        i += 1
    end
    return groups
end

function parse_group(lines::Vector{String})
    text = join(lines, " ")
    tokens = tokenize(text)
    # Remove any stray '&' tokens (they were group delimiters)
    filter!(t -> t != "&", tokens)

    dict = Dict{String, Any}()
    eq_positions = findall(t -> t == "=", tokens)
    isempty(eq_positions) && return dict

    for (idx, eq_idx) in enumerate(eq_positions)
        # ---- LHS ----
        lhs_tokens = String[]
        j = eq_idx - 1
        while j >= 1 && (is_variable_name(tokens[j]) || tokens[j] == "(" ||
                         tokens[j] == ")" || (tryparse(Int, tokens[j]) !== nothing))
            pushfirst!(lhs_tokens, tokens[j])
            j -= 1
        end
        if isempty(lhs_tokens)
            error("Empty LHS before '=' at position $eq_idx")
        end
        base_name, index = parse_lhs(lhs_tokens)

        # ---- RHS ----
        next_eq = (idx < length(eq_positions)) ? eq_positions[idx+1] : length(tokens) + 1
        rhs_tokens = tokens[eq_idx+1 : next_eq-1]
        value = parse_value(rhs_tokens)

        # Store
        if index === nothing
            dict[base_name] = value
        else
            if !haskey(dict, base_name) || !(dict[base_name] isa Dict)
                dict[base_name] = Dict{Int,Any}()
            end
            dict[base_name][index] = value
        end
    end

    return dict
end

function parse_lhs(tokens::Vector{String})
    if length(tokens) == 1
        return tokens[1], nothing
    elseif length(tokens) >= 4 && tokens[2] == "(" && tokens[end] == ")"
        idx_token = tokens[3]
        if tryparse(Int, idx_token) !== nothing
            return tokens[1], parse(Int, idx_token)
        else
            error("Invalid array index: $idx_token")
        end
    else
        error("Unsupported LHS expression: $(join(tokens, ""))")
    end
end

function tokenize(text::String)
    tokens = Vector{String}()
    i = 1
    n = length(text)
    while i <= n
        c = text[i]
        if isspace(c)
            i += 1
            continue
        elseif c == '\'' || c == '"'
            quote_char = c
            i += 1
            start = i
            while i <= n && text[i] != quote_char
                i += 1
            end
            if i > n
                error("Unclosed string literal")
            end
            str = text[start:i-1]
            push!(tokens, str)
            i += 1
        elseif c == '=' || c == ',' || c == '(' || c == ')' || c == '&'
            push!(tokens, string(c))
            i += 1
        else
            start = i
            while i <= n && !isspace(text[i]) && text[i] != '=' &&
                  text[i] != ',' && text[i] != '(' && text[i] != ')' && text[i] != '&'
                i += 1
            end
            token = text[start:i-1]
            push!(tokens, token)
        end
    end
    return tokens
end

function parse_value(tokens::Vector{String})
    if any(t -> t == ",", tokens)
        values = []
        current = String[]
        for tok in tokens
            if tok == ","
                if !isempty(current)
                    push!(values, parse_scalar(join(current, "")))
                    current = []
                end
            else
                push!(current, tok)
            end
        end
        if !isempty(current)
            push!(values, parse_scalar(join(current, "")))
        end
        return values
    else
        return parse_scalar(join(tokens, ""))
    end
end

function parse_scalar(str::String)
    s = strip(str)
    isempty(s) && return nothing
    # Booleans
    if s in (".true.", ".TRUE.", "true", "TRUE", "T", "t")
        return true
    elseif s in (".false.", ".FALSE.", "false", "FALSE", "F", "f")
        return false
    end
    # Integer
    if tryparse(Int, s) !== nothing
        return parse(Int, s)
    end
    # Float
    if tryparse(Float64, s) !== nothing
        return parse(Float64, s)
    end
    # Otherwise treat as string (unquoted)
    return s
end

is_variable_name(s::String) = occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", s)