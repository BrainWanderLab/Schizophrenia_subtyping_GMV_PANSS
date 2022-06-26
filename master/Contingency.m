    function Cont=Contingency(Mem1,Mem2)
                      
        Cont=zeros(max(Mem1),max(Mem2));
        
        for i = 1:length(Mem1);
            Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
        end
    end