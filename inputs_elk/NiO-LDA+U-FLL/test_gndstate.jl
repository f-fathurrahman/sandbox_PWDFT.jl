using ElkDFTWrapper
const elk = ElkDFTWrapper

#elk.init_run() # call read_input, init0, init1
elk.call_read_input()
elk.call_gndstate() # call the original gndstate

