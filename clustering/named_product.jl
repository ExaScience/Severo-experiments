module NamedProduct
    struct NamedProductIterator{names,T<:Tuple}
        iterators::NamedTuple{names, T}
    end

    named_product(;kv...) = NamedProductIterator(values(kv))

    import Base: IteratorSize, size, axes, length, HasEltype, HasShape, IsInfinite, HasLength, tail
    IteratorSize(::Type{NamedProductIterator{(), Tuple{}}}) = HasShape{0}()
    function IteratorSize(::Type{NamedProductIterator{names, T}}) where {names, T<:Tuple}
        prod_iteratorsize(
            IteratorSize(Base.tuple_type_head(T)),
            IteratorSize(NamedProductIterator{tail(names),Base.tuple_type_tail(T)}) )
    end

    prod_iteratorsize(::HasLength, ::HasLength) = HasShape{2}()
    prod_iteratorsize(::HasLength, ::HasShape{N}) where {N} = HasShape{N+1}()
    prod_iteratorsize(::HasShape{N}, ::HasLength) where {N} = HasShape{N+1}()
    prod_iteratorsize(::HasShape{M}, ::HasShape{N}) where {M,N} = HasShape{M+N}()

    # products can have an infinite iterator
    prod_iteratorsize(::IsInfinite, ::IsInfinite) = IsInfinite()
    prod_iteratorsize(a, ::IsInfinite) = IsInfinite()
    prod_iteratorsize(::IsInfinite, b) = IsInfinite()
    prod_iteratorsize(a, b) = SizeUnknown()

    axes(P::NamedProductIterator) = _prod_indices(P.iterators)
    _prod_indices(::NamedTuple{(), Tuple{}}) = ()
    _prod_indices(t::NamedTuple) = (_prod_axes1(t[1], IteratorSize(t[1]))..., _prod_indices(tail(t))...)
    _prod_axes1(a, ::HasShape)  = axes(a)
    _prod_axes1(a, ::HasLength) = (OneTo(length(a)),)
    _prod_axes1(a, A) =
        throw(ArgumentError("Cannot compute indices for object of type $(typeof(a))"))

    size(P::NamedProductIterator) = _prod_size(P.iterators)
    _prod_size(::NamedTuple{(), Tuple{}}) = ()
    _prod_size(t::NamedTuple) = (_prod_size1(t[1], IteratorSize(t[1]))..., _prod_size(tail(t))...)
    _prod_size1(a, ::HasShape)  = size(a)
    _prod_size1(a, ::HasLength) = (length(a),)
    _prod_size1(a, A) =
        throw(ArgumentError("Cannot compute size for object of type $(typeof(a))"))

    length(P::NamedProductIterator) = prod(size(P))

    import Base: IteratorEltype, eltype, EltypeUnknown, iterate, isdone
    IteratorEltype(::Type{NamedProductIterator{(),Tuple{}}}) = HasEltype()
    IteratorEltype(::Type{NamedProductIterator{Tuple{I}}}) where {I} = IteratorEltype(I)
    function IteratorEltype(::Type{NamedProductIterator{names, T}}) where {names, T<:Tuple}
        I = Base.tuple_type_head(T)
        P = NamedProductIterator{Base.tail(names), Base.tuple_type_tail(T)}
        IteratorEltype(I) == EltypeUnknown() ? EltypeUnknown() : IteratorEltype(P)
    end

    eltype(::Type{NamedProductIterator{names, I}}) where {names, I} = NamedTuple{names, _prod_eltype(I)}
    _prod_eltype(::Type{Tuple{}}) = Tuple{}
    _prod_eltype(::Type{I}) where {I<:Tuple} =
        Base.tuple_type_cons(eltype(Base.tuple_type_head(I)),_prod_eltype(Base.tuple_type_tail(I)))

    iterate(::NamedProductIterator{(),Tuple{}}) = (), true
    iterate(::NamedProductIterator{(),Tuple{}}, state) = nothing

    @inline isdone(P::NamedProductIterator) = any(isdone, P.iterators)
    @inline function _pisdone(iters, states)
        iter1 = first(iters)
        done1 = isdone(iter1, first(states)[2]) # check step
        done1 === true || return done1 # false or missing
        done1 = isdone(iter1) # check restart
        done1 === true || return done1 # false or missing
        return _pisdone(tail(iters), tail(states)) # check tail
    end
    @inline isdone(P::NamedProductIterator, states) = _pisdone(P.iterators, states)

    @inline _piterate() = ()
    @inline function _piterate(iter1, rest...)
        next = iterate(iter1)
        next === nothing && return nothing
        restnext = _piterate(rest...)
        restnext === nothing && return nothing
        return (next, restnext...)
    end

    @inline function iterate(P::NamedProductIterator{names,T}) where {names, T}
        isdone(P) === true && return nothing
        next = _piterate(P.iterators...)
        next === nothing && return nothing
        return (NamedTuple{names}(map(first, next)), next)
    end

    @inline _piterate1(::NamedTuple{(),Tuple{}}, ::Tuple{}) = nothing
    @inline function _piterate1(iters, states)
        iter1 = first(iters)
        next = iterate(iter1, first(states)[2])
        restnext = tail(states)
        if next === nothing
            isdone(iter1) === true && return nothing
            restnext = _piterate1(tail(iters), restnext)
            restnext === nothing && return nothing
            next = iterate(iter1)
            next === nothing && return nothing
        end
        return (next, restnext...)
    end

    @inline function iterate(P::NamedProductIterator{names,T}, states) where {names,T}
        isdone(P, states) === true && return nothing
        next = _piterate1(P.iterators, states)
        next === nothing && return nothing
        return (NamedTuple{names}(map(first, next)), next)
    end

    export NamedProductIterator, named_product
end
