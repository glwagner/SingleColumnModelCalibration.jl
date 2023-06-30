sleepy() = sleep(1)

start = time_ns()

#@sync begin
#@sync for i = 1:10
for i = 1:10
        @async sleepy()
        #sleepy()
    end
#end

stop = time_ns()

@info 1e-9 * (stop - start)


