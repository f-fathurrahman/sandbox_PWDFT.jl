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
    dict = Dict{String, Any}()
    i = 1
    while i <= length(tokens)
        tok = tokens[i]
        if tok isa String && tok != "=" && tok != "(" && tok != ")" && tok != "," && tok != "/" && tok != "&"
            name = tok
            index = nothing
            # check for array element
            if i+1 <= length(tokens) && tokens[i+1] == "("
                i += 2
                if i > length(tokens)
                    error("Expected index after '('")
                end
                idx_token = tokens[i]
                if tryparse(Int, idx_token) === nothing
                    error("Expected integer index, got '$idx_token'")
                end
                index = parse(Int, idx_token)
                i += 1
                if i > length(tokens) || tokens[i] != ")"
                    error("Expected ')' after index")
                end
                i += 1
                if i > length(tokens) || tokens[i] != "="
                    error("Expected '=' after ')'")
                end
                i += 1
                val_tokens = String[]
                while i <= length(tokens)
                    t = tokens[i]
                    if is_variable_name(t) && !is_boolean_literal(t) && t != name
                        break
                    else
                        push!(val_tokens, t)
                        i += 1
                    end
                end
                value = parse_value(val_tokens)
                if !haskey(dict, name) || !(dict[name] isa Dict)
                    dict[name] = Dict{Int,Any}()
                end
                dict[name][index] = value
            else
                # scalar assignment
                i += 1
                if i > length(tokens) || tokens[i] != "="
                    println("tokens = $tokens")
                    error("Expected '=' after variable name '$name'")
                end
                i += 1
                val_tokens = String[]
                while i <= length(tokens)
                    t = tokens[i]
                    if is_variable_name(t) && !is_boolean_literal(t)
                        break
                    else
                        push!(val_tokens, t)
                        i += 1
                    end
                end
                value = parse_value(val_tokens)
                dict[name] = value
            end
        else
            i += 1
        end
    end
    return dict
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
        elseif c == '=' || c == ',' || c == '/' || c == '&' || c == '(' || c == ')'
            push!(tokens, string(c))
            i += 1
        elseif c == '\'' || c == '"'
            quote_char = c
            i += 1
            start = i
            while i <= n && text[i] != quote_char
                i += 1
            end
            if i > n
                println("text = $text c = $c")
                error("Unclosed string literal")
            end
            str = text[start:i-1]
            push!(tokens, str)
            i += 1
        elseif c == '.' || c == '-' || c == '+' || isdigit(c)
            start = i
            while i <= n && (isdigit(text[i]) || text[i] == '.' ||
                             text[i] == 'e' || text[i] == 'E' ||
                             text[i] == '-' || text[i] == '+')
                i += 1
            end
            token = text[start:i-1]
            push!(tokens, token)
        elseif isletter(c) || c == '_'
            start = i
            while i <= n && (isletter(text[i]) || isdigit(text[i]) || text[i] == '_')
                i += 1
            end
            token = text[start:i-1]
            push!(tokens, token)
        else
            # skip any other unknown characters
            i += 1
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
    if isempty(s)
        return nothing
    end
    if s in (".true.", ".TRUE.", "true", "TRUE", "T", "t")
        return true
    elseif s in (".false.", ".FALSE.", "false", "FALSE", "F", "f")
        return false
    end
    if tryparse(Int, s) !== nothing
        return parse(Int, s)
    end
    if tryparse(Float64, s) !== nothing
        return parse(Float64, s)
    end
    return s
end

is_variable_name(s::String) = occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", s)
is_boolean_literal(s::String) = s in (".true.", ".TRUE.", "true", "TRUE", "T", "t",
                                      ".false.", ".FALSE.", "false", "FALSE", "F", "f")