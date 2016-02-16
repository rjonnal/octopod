
            prof = np.mean(im2,axis=1)
            newprof = np.zeros(prof.shape)
            offset = (-1)*ty
            print offset
            if offset<0:
                newprof[-offset:] = prof[:offset]
                if False:
                    plt.figure(2)
                    plt.plot(self.prof/count)
                    plt.plot(newprof)
                    plt.pause(.5)
            elif offset>0:
                newprof[:-offset] = prof[offset:]
                if False:
                    plt.figure(2)
                    plt.plot(self.prof/count)
                    plt.plot(newprof)
                    plt.show(.5)
            else:
                newprof = prof
            self.prof = self.prof + newprof
            count = count + 1
            
            if debug:
                plt.figure(1)
                plt.cla()
                plt.plot(self.prof)
                #plt.imshow(im1,aspect='auto',interpolation='none')
                plt.pause(.0001)
