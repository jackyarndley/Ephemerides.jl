
export ephem_vector3, ephem_vector6, ephem_vector9, ephem_vector12, 
       ephem_rotation3, ephem_rotation6, ephem_rotation9, ephem_rotation12

for (order, pfun1, afun1, pfun2, afun2) in zip(
    (1, 2, 3, 4),
    (:ephem_vector3, :ephem_vector6, :ephem_vector9, :ephem_vector12),
    (:ephem_rotation3, :ephem_rotation6, :ephem_rotation9, :ephem_rotation12),
    (:spk_vector3, :spk_vector6, :spk_vector9, :spk_vector12),
    (:pck_vector3, :pck_vector6, :pck_vector9, :pck_vector12)
)

    @eval begin 

        """
            $($pfun1)(eph::EphemerisProvider, from::Int, to::Int, time::Number)

        Compute the $(3*$order)-elements state of one body (to) relative to another (from)
        at `time`, expressed in TDB/TCB seconds since J2000, in accordance with the kernel
        timescale.
        """
        function ($pfun1)(eph::EphemerisProvider, from::Int, to::Int, time::Number)

            links = spk_links(eph)
            if haskey(links, to)
                to_links = links[to]
                if haskey(to_links, from)
                    files = get_daf(eph)
                    for link in to_links[from] 
                        if initial_time(link) <= time <= final_time(link)   
                            return factor(link)*$(pfun2)(files[file_id(link)], link, time)
                        end
                    end
                end
            else 
                throw(
                    jEph.EphemerisError(
                        "ephemeris data for point with NAIFId $to with respect to point " * 
                        "$from is unavailable."
                    )
                )
            end
        
            if !haskey(links[to], from)
                throw(
                    jEph.EphemerisError(
                        "ephemeris data for point with NAIFId $to with respect to point " * 
                        "$from is unavailable."
                    )
                )
            end

            throw(
                jEph.EphemerisError(
                    "ephemeris data for point with NAIFId $to with respect to point " *
                    "$from is not available at $(time) seconds since J2000."
                ),
            )

        end

        """
            $($afun1)(eph::EphemerisProvider, from::Int, to::Int, time::Number)

        Compute the $(3*$order)-elements orientation angles of one set of axes (to) relative 
        to another (from) at `time`, expressed in TDB/TCB seconds since J2000, in accordance 
        with the kernel timescale.

        !!! note 
            For the orientation angles, it is not possible to automatically compute the 
            reverse transformation , i.e., if the orientation of PA440 is defined 
            with respect to the ICRF, it is not possible to compute the rotation from the 
            PA440 to the ICRF with this routine.
        """
        function ($afun1)(eph::EphemerisProvider, from::Int, to::Int, time::Number)
            links = pck_links(eph)
            if haskey(links, to)
                to_links = links[to]
                if haskey(to_links, from)
                    files = get_daf(eph)
                    for link in to_links[from] 
                        if initial_time(link) <= time <= final_time(link)   
                            return $(afun2)(files[file_id(link)], link, time)
                        end
                    end
                end
            else 
                throw(
                    jEph.EphemerisError(
                        "ephemeris data for axes with NAIFId $to with respect to axes " * 
                        "$from is unavailable."
                    )
                )
            end
        
            if !haskey(links[to], from)
                throw(
                    jEph.EphemerisError(
                        "ephemeris data for axes with NAIFId $to with respect to axes " * 
                        "$from is unavailable."
                    )
                )
            end
        
            throw(
                jEph.EphemerisError(
                    "ephemeris data for axes with NAIFId $to with respect to axes " *
                    "$from is not available at $(time) seconds since J2000."
                ),
            )

        end

        # This is internal for the SPKs
        function ($pfun2)(daf::DAF, link::SPKLink, time::Number)

            # Retrieve list and element link IDs
            lid = link.lid
            eid = link.eid
            seglist = segment_list(daf)
            
            if lid == 1 
                return $(pfun2)(daf, seglist.spk2[eid], time)
            elseif lid == 2
                return $(pfun2)(daf, seglist.spk9[eid], time)
            elseif lid == 3
                return $(pfun2)(daf, seglist.spk1[eid], time)
            elseif lid == 4
                return $(pfun2)(daf, seglist.spk14[eid], time)
            elseif lid == 5
                return $(pfun2)(daf, seglist.spk15[eid], time)
            elseif lid == 6
                return $(pfun2)(daf, seglist.spk8[eid], time)
            elseif lid == 7
                return $(pfun2)(daf, seglist.spk19[eid], time)
            elseif lid == 8
                return $(pfun2)(daf, seglist.spk20[eid], time)
            elseif lid == 9
                return $(pfun2)(daf, seglist.spk5[eid], time)
            else
                return $(pfun2)(daf, seglist.spk17[eid], time)
            end
        end

        # This is internal for the PCKs
        function ($afun2)(daf::DAF, link::SPKLink, time::Number)

            # Retrieve list and element link IDs
            lid = link.lid
            eid = link.eid
            seglist = segment_list(daf)
            
            if lid == 1 
                return $(pfun2)(daf, seglist.spk2[eid], time)
            else
                return $(pfun2)(daf, seglist.spk20[eid], time)
            end
 
        end

    end

end
