require 'listen'
require 'terminal-notifier'

TerminalNotifier.notify('Hello World', :open => 'http://apple.com')


Listen.to!('.', :relative_paths => false ) do |modified, added, removed|
	if (!modified.empty?)  
		puts  modified.inspect
	end
	 
	if (!added.empty?) 
		puts  added[0]
		TerminalNotifier.notify('Hello World',  :title => 'Ruby', :subtitle => 'Programming Language', :open => "file://#{added[0]}")

	end

 	if (!removed.empty?) 
 		puts  removed.inspect
 	end
end


